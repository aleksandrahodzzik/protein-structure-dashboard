"""Parse mmCIF into protein chains: observed sequence, residues, ligands and SEQRES."""

from __future__ import annotations

import io
from dataclasses import dataclass
from itertools import groupby
from typing import Any

from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1_extended

from .analysis import STANDARD_AMINO_ACIDS

UNKNOWN_RESIDUE = "X"
WATER_NAMES = frozenset({"HOH", "WAT", "DOD"})
BACKBONE_ATOMS = frozenset({"N", "CA", "C"})
# What a 3Dmol cartoon needs: the protein backbone (CA alone draws nothing) and the
# sugar-phosphate backbone of nucleic acids.
CARTOON_ATOMS = frozenset({"N", "CA", "C", "O", "P", "O5'", "C5'", "C4'", "C3'", "O3'"})
# Categories kept in backbone-only mode so helices and sheets are still drawn as such.
SECONDARY_STRUCTURE_CATEGORIES = frozenset({"_struct_conf", "_struct_sheet_range"})
# Two residues are neighbours in the chain if the peptide bond C(i)-N(i+1) is present
# (1.33 A; 1.8 A is Biopython PPBuilder's cut-off). Models with only CA atoms have no
# bond to measure, so CA(i)-CA(i+1) (3.8 A; 4.3 A is CaPPBuilder's cut-off) is used.
PEPTIDE_BOND_MAX = 1.8
CA_CA_MAX = 4.3
# The deposited sequence and how each position maps onto the author numbering.
SEQRES_CATEGORY = "_pdbx_poly_seq_scheme"
NO_INSERTION_CODE = frozenset({".", "?"})
RANGE_DASH = "\u2013"  # en dash, the typeset range sign


class StructureParseError(ValueError):
    """The text is not usable mmCIF, or it contains no protein."""


@dataclass(frozen=True)
class ResidueRecord:
    number: int
    insertion_code: str
    name: str  # three-letter code as written in the file, e.g. "MSE"
    letter: str  # one-letter code of the parent residue, "X" if unknown
    segment: int = 0  # index of the contiguous modelled stretch the residue belongs to

    @property
    def label(self) -> str:
        return f"{self.number}{self.insertion_code}"


@dataclass(frozen=True)
class SeqresResidue:
    """One position of the deposited (SEQRES) sequence, in author numbering."""

    number: int
    insertion_code: str
    names: tuple[str, ...]  # usually one; several where the deposited sequence is heterogeneous

    @property
    def name(self) -> str:
        return self.names[0]

    @property
    def letter(self) -> str:
        return residue_letter(self.name)

    @property
    def label(self) -> str:
        return f"{self.number}{self.insertion_code}"


@dataclass(frozen=True)
class Ligand:
    name: str
    number: int
    chain_id: str


@dataclass(frozen=True)
class ProteinChain:
    chain_id: str
    residues: tuple[ResidueRecord, ...]
    ligands: tuple[Ligand, ...]
    seqres: tuple[SeqresResidue, ...] = ()  # empty when the file has no SEQRES mapping

    @property
    def sequence(self) -> str:
        """Observed (modelled) sequence; may contain X for unknown residues."""
        return "".join(r.letter for r in self.residues)

    @property
    def residue_count(self) -> int:
        return len(self.residues)

    @property
    def unknown_count(self) -> int:
        return sum(r.letter == UNKNOWN_RESIDUE for r in self.residues)

    @property
    def segment_count(self) -> int:
        """Contiguous modelled stretches; more than one means gaps in the density."""
        return len({r.segment for r in self.residues})

    @property
    def breaks(self) -> tuple[int, ...]:
        """Indices ``i`` where residue ``i`` and ``i + 1`` are not bonded neighbours."""
        residues = self.residues
        return tuple(
            i for i in range(len(residues) - 1) if residues[i].segment != residues[i + 1].segment
        )

    @property
    def ligand_names(self) -> tuple[str, ...]:
        return tuple(sorted({lig.name for lig in self.ligands}))


@dataclass(frozen=True)
class MissingStretch:
    """Consecutive SEQRES positions with no modelled residue."""

    first: SeqresResidue
    last: SeqresResidue
    count: int

    @property
    def label(self) -> str:
        if self.count == 1:
            return self.first.label
        # A minus sign next to the range dash ("-18" + dash + "26") is hard to read, so a
        # range that starts below zero is written with words instead.
        joiner = " to " if self.first.number < 0 else RANGE_DASH
        return f"{self.first.label}{joiner}{self.last.label}"


@dataclass(frozen=True)
class Mismatch:
    """A modelled residue whose type differs from the SEQRES residue at its position."""

    label: str
    seqres_name: str
    modelled_name: str


@dataclass(frozen=True)
class SeqresComparison:
    """Where the observed sequence of a chain differs from its deposited sequence."""

    seqres: tuple[SeqresResidue, ...]
    modelled: tuple[ResidueRecord | None, ...]  # per SEQRES position, None if not modelled
    missing: tuple[MissingStretch, ...]
    mismatches: tuple[Mismatch, ...]
    unmatched: tuple[ResidueRecord, ...]  # modelled residues with no SEQRES position

    @property
    def seqres_length(self) -> int:
        return len(self.seqres)

    @property
    def modelled_count(self) -> int:
        return sum(r is not None for r in self.modelled)

    @property
    def missing_count(self) -> int:
        return self.seqres_length - self.modelled_count

    def alignment(self, width: int = 60, block: int = 10) -> str:
        """SEQRES above the modelled residues, ``-`` where a residue was not modelled."""
        starts = range(0, self.seqres_length, width)
        label_width = max(len(self.seqres[start].label) for start in starts)
        rows: list[str] = []
        for start in starts:
            top = "".join(s.letter for s in self.seqres[start : start + width])
            bottom = "".join(
                "-" if r is None else r.letter for r in self.modelled[start : start + width]
            )
            label = self.seqres[start].label
            rows.append(
                f"SEQRES {label:>{label_width}} {_blocks(top, block)}\n"
                f"Model  {'':>{label_width}} {_blocks(bottom, block)}"
            )
        return "\n\n".join(rows)


def _blocks(letters: str, size: int) -> str:
    return " ".join(letters[i : i + size] for i in range(0, len(letters), size))


def compare_to_seqres(chain: ProteinChain) -> SeqresComparison | None:
    """Line the modelled residues up with the deposited sequence by author numbering.

    Returns None when the file has no SEQRES mapping for this chain.
    """
    if not chain.seqres:
        return None
    observed = {(r.number, r.insertion_code): r for r in chain.residues}
    modelled = tuple(observed.get((s.number, s.insertion_code)) for s in chain.seqres)
    missing = []
    position = 0
    for is_missing, run in groupby(modelled, key=lambda r: r is None):
        length = len(list(run))
        if is_missing:
            first, last = chain.seqres[position], chain.seqres[position + length - 1]
            missing.append(MissingStretch(first=first, last=last, count=length))
        position += length
    mismatches = tuple(
        Mismatch(label=s.label, seqres_name=s.name, modelled_name=r.name)
        for s, r in zip(chain.seqres, modelled, strict=True)
        if r is not None and r.name not in s.names
    )
    in_seqres = {(s.number, s.insertion_code) for s in chain.seqres}
    return SeqresComparison(
        seqres=chain.seqres,
        modelled=modelled,
        missing=tuple(missing),
        mismatches=mismatches,
        unmatched=tuple(r for r in chain.residues if (r.number, r.insertion_code) not in in_seqres),
    )


@dataclass(frozen=True)
class ParsedStructure:
    pdb_id: str
    chains: tuple[ProteinChain, ...]
    atom_count: int

    def chain(self, chain_id: str) -> ProteinChain:
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain
        raise KeyError(chain_id)


def residue_letter(resname: str) -> str:
    """One-letter code for a residue name.

    Modified residues map to their parent amino acid (MSE selenomethionine -> M,
    SEP phosphoserine -> S, ...) using Biopython's table; anything it does not
    know (including UNK) becomes X rather than a guess. Letters outside the 20
    standard amino acids also become X, so the sequence and the analysis module
    always agree on which residues are usable.
    """
    letter = protein_letters_3to1_extended.get(resname.strip().upper(), UNKNOWN_RESIDUE)
    return letter if letter in STANDARD_AMINO_ACIDS else UNKNOWN_RESIDUE


def _is_polymer_residue(residue: Any) -> bool:
    """True for amino-acid residues, including modified ones written as HETATM.

    Standard names always count. Anything else needs a full N/CA/C backbone, which
    keeps ions such as calcium (residue and atom both named CA) out of the sequence.
    """
    if residue.get_resname().strip() in WATER_NAMES:
        return False
    if is_aa(residue, standard=True):
        return True
    return residue.child_dict.keys() >= BACKBONE_ATOMS


def _atoms(residue: Any, name: str) -> list[Any]:
    """Every alternative position of an atom (one, unless it is disordered)."""
    atom = residue[name]
    return list(atom.disordered_get_list()) if atom.is_disordered() else [atom]


def _within(a: Any, b: Any, name_a: str, name_b: str, cutoff: float) -> bool:
    return any(x - y < cutoff for x in _atoms(a, name_a) for y in _atoms(b, name_b))


def _bonded(previous: Any, current: Any) -> bool:
    """True if two polymer residues are neighbours in the chain, not across a gap."""
    if "C" in previous and "N" in current:
        return _within(previous, current, "C", "N", PEPTIDE_BOND_MAX)
    if "CA" in previous and "CA" in current:
        return _within(previous, current, "CA", "CA", CA_CA_MAX)
    return False


def _protein_residues(chain: Any) -> list[tuple[Any, int]]:
    """Polymer residues in file order, each with the index of its contiguous segment.

    A new segment starts wherever a residue is not bonded to the previous one, which
    is what residues missing from the electron density look like in the coordinates.
    A lone residue between two disordered loops is a segment of its own.
    """
    residues = [res for res in chain if _is_polymer_residue(res)]
    segments: list[tuple[Any, int]] = []
    segment = -1
    for index, res in enumerate(residues):
        if index == 0 or not _bonded(residues[index - 1], res):
            segment += 1
        segments.append((res, segment))
    return segments


def _category_text(mmcif_text: str, category: str) -> str | None:
    """The lines of one mmCIF category, found by text search instead of parsing the
    whole (possibly very large) file. RCSB ends every category with a '#' line."""
    start = mmcif_text.find(f"\n{category}.")
    if start < 0:
        return None
    loop = mmcif_text.rfind("\nloop_", 0, start)
    if loop >= 0 and mmcif_text[loop:start].strip() == "loop_":
        start = loop
    end = mmcif_text.find("\n#", start + 1)
    return mmcif_text[start + 1 : end if end >= 0 else len(mmcif_text)]


def _seqres_by_chain(mmcif_text: str) -> dict[str, tuple[SeqresResidue, ...]]:
    """Deposited sequence of every polymer chain, keyed by author chain ID.

    ``_pdbx_poly_seq_scheme`` lists every SEQRES position with its author number, also
    for residues that have no coordinates, which is what lets missing residues be named.
    A position listed twice is a heterogeneous one (two residue types deposited).
    """
    text = _category_text(mmcif_text, SEQRES_CATEGORY)
    if text is None:
        return {}
    try:
        table = MMCIF2Dict(io.StringIO(f"data_seqres\n{text}"))
        columns = [
            table[f"{SEQRES_CATEGORY}.{name}"]
            for name in ("pdb_strand_id", "asym_id", "seq_id", "mon_id", "pdb_seq_num")
        ]
        insertion_codes = table.get(f"{SEQRES_CATEGORY}.pdb_ins_code", ["."] * len(columns[0]))
    except (KeyError, ValueError):  # malformed or incomplete: carry on without SEQRES
        return {}
    if len({len(column) for column in [*columns, insertion_codes]}) != 1:
        return {}  # a broken loop with ragged columns: no SEQRES rather than no structure

    positions: dict[str, dict[tuple[str, str], tuple[int, str, list[str]]]] = {}
    for (strand, asym, seq_id, name, number), code in zip(
        zip(*columns, strict=True), insertion_codes, strict=True
    ):
        if not number.lstrip("-").isdigit():
            continue
        chain = positions.setdefault(strand, {})
        key = (asym, seq_id)
        if key not in chain:
            chain[key] = (int(number), "" if code in NO_INSERTION_CODE else code, [name])
        elif name not in chain[key][2]:
            chain[key][2].append(name)
    return {
        strand: tuple(
            SeqresResidue(number=number, insertion_code=code, names=tuple(names))
            for number, code, names in chain.values()
        )
        for strand, chain in positions.items()
    }


def trim_mmcif(mmcif_text: str, *, backbone_only: bool = False) -> str:
    """Rewrite the mmCIF so the viewer draws what the analysis analysed.

    Only the first model is kept (an NMR ensemble would otherwise be drawn as ten
    overlapping copies). With ``backbone_only`` just the atoms a cartoon needs and the
    helix/sheet records are kept, which shrinks huge entries by an order of magnitude.
    Works on the text instead of re-parsing, which would double the cost on big files.
    """
    blocks = _category_blocks(mmcif_text)
    trimmed: list[list[str]] = []
    found_atoms = False
    for block in blocks:
        category = _block_category(block)
        if category == "_atom_site":
            found_atoms = True
            block = _filter_atom_rows(block, backbone_only)
        elif backbone_only and category not in SECONDARY_STRUCTURE_CATEGORIES | {None}:
            continue
        trimmed.append(block)
    if not found_atoms:
        raise StructureParseError("mmCIF data has no atom_site loop.")
    return "\n#\n".join("\n".join(block) for block in trimmed) + "\n#\n"


def _category_blocks(mmcif_text: str) -> list[list[str]]:
    """Split on the '#' separator lines RCSB writes between mmCIF categories."""
    blocks: list[list[str]] = [[]]
    for line in mmcif_text.splitlines():
        if line.startswith("#"):
            blocks.append([])
        else:
            blocks[-1].append(line)
    return [block for block in blocks if block]


def _block_category(block: list[str]) -> str | None:
    """'_atom_site' for a block of _atom_site.* items; None for the data_ header."""
    for line in block:
        if line.startswith("_"):
            return line.split(".")[0].split()[0]
    return None


def _filter_atom_rows(block: list[str], backbone_only: bool) -> list[str]:
    first_item = next(i for i, line in enumerate(block) if line.startswith("_"))
    header_end = next(i for i in range(first_item, len(block)) if not block[i].startswith("_"))
    columns = [line.split()[0] for line in block[first_item:header_end]]
    name_column = columns.index("_atom_site.label_atom_id")
    residue_column = columns.index("_atom_site.label_comp_id")
    model_column = columns.index("_atom_site.pdbx_PDB_model_num")

    kept = block[:header_end]
    first_model: str | None = None
    for row in block[header_end:]:
        tokens = row.split()
        first_model = first_model or tokens[model_column]
        if tokens[model_column] != first_model:
            continue
        if backbone_only and (
            tokens[name_column].strip("\"'") not in CARTOON_ATOMS
            or tokens[residue_column] in WATER_NAMES  # water oxygens are also named "O"
        ):
            continue
        kept.append(row)
    return kept


def parse_structure(mmcif_text: str, pdb_id: str = "") -> ParsedStructure:
    """Parse the first model of an mmCIF file into protein chains.

    Chain IDs are the author chain IDs, which is what RCSB shows on entry pages.
    Chains without amino acids (DNA/RNA, ligand-only) are omitted.
    """
    try:
        structure = MMCIFParser(QUIET=True).get_structure(
            pdb_id or "structure", io.StringIO(mmcif_text)
        )
        model = structure[0]
    except Exception as exc:  # Biopython raises a mix of KeyError/ValueError/IndexError
        raise StructureParseError(f"Could not parse mmCIF data: {exc}") from exc

    seqres = _seqres_by_chain(mmcif_text)
    chains: list[ProteinChain] = []
    for chain in model:
        residues = _protein_residues(chain)
        if not residues:
            continue
        polymer_ids = {res.id for res, _ in residues}
        ligands = tuple(
            Ligand(name=res.get_resname().strip(), number=res.id[1], chain_id=chain.id)
            for res in chain
            if res.id[0] not in (" ", "W")
            and res.id not in polymer_ids
            and res.get_resname().strip() not in WATER_NAMES
        )
        chains.append(
            ProteinChain(
                chain_id=chain.id,
                residues=tuple(
                    ResidueRecord(
                        number=res.id[1],
                        insertion_code=res.id[2].strip(),
                        name=res.get_resname().strip(),
                        letter=residue_letter(res.get_resname()),
                        segment=segment,
                    )
                    for res, segment in residues
                ),
                ligands=ligands,
                seqres=seqres.get(chain.id, ()),
            )
        )

    if not chains:
        raise StructureParseError("This entry contains no protein chains to analyse.")
    return ParsedStructure(
        pdb_id=pdb_id,
        chains=tuple(chains),
        atom_count=sum(1 for _ in model.get_atoms()),
    )
