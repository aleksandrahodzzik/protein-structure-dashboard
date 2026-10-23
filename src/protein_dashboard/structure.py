"""Parse mmCIF into protein chains: observed sequence, residues and ligands."""

from __future__ import annotations

import io
import warnings
from dataclasses import dataclass
from typing import Any

from Bio.PDB import MMCIFParser, PPBuilder
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


class StructureParseError(ValueError):
    """The text is not usable mmCIF, or it contains no protein."""


@dataclass(frozen=True)
class ResidueRecord:
    number: int
    insertion_code: str
    name: str  # three-letter code as written in the file, e.g. "MSE"
    letter: str  # one-letter code of the parent residue, "X" if unknown


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
    segment_count: int  # contiguous modelled stretches; > 1 means gaps in the density

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
    def ligand_names(self) -> tuple[str, ...]:
        return tuple(sorted({lig.name for lig in self.ligands}))


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


def _protein_residues(chain: Any) -> tuple[list[Any], int]:
    """Polymer residues in file order, and how many contiguous segments they form.

    PPBuilder links residues by C-N distance, which gives us chain breaks, but it
    only emits runs of two or more. Isolated residues (e.g. a lone residue between
    two disordered loops) are therefore added back as segments of their own so the
    sequence stays complete. aa_only=0 lets it link modified residues like MSE;
    its "unknown modified amino acid" warning fires for every ion and is silenced
    because _is_polymer_residue does the filtering explicitly.
    """
    residues = [res for res in chain if _is_polymer_residue(res)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        peptides = PPBuilder().build_peptides(chain, aa_only=0)
    linked = {res.id for peptide in peptides for res in peptide}
    runs = sum(_is_polymer_residue(peptide[0]) for peptide in peptides)
    isolated = sum(res.id not in linked for res in residues)
    return residues, runs + isolated


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

    chains: list[ProteinChain] = []
    for chain in model:
        residues, segments = _protein_residues(chain)
        if not residues:
            continue
        polymer_ids = {res.id for res in residues}
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
                    )
                    for res in residues
                ),
                ligands=ligands,
                segment_count=segments,
            )
        )

    if not chains:
        raise StructureParseError("This entry contains no protein chains to analyse.")
    return ParsedStructure(
        pdb_id=pdb_id,
        chains=tuple(chains),
        atom_count=sum(1 for _ in model.get_atoms()),
    )
