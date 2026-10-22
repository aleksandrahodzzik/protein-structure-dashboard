"""Parse mmCIF into protein chains: observed sequence, residues and ligands."""

from __future__ import annotations

import io
import warnings
from dataclasses import dataclass
from typing import Any

from Bio.PDB import MMCIFParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1_extended

UNKNOWN_RESIDUE = "X"
WATER_NAMES = frozenset({"HOH", "WAT", "DOD"})
BACKBONE_ATOMS = frozenset({"N", "CA", "C"})
BACKBONE_TRACE_ATOMS = frozenset({"CA", "P"})


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
    know (including UNK) becomes X rather than a guess.
    """
    return protein_letters_3to1_extended.get(resname.strip().upper(), UNKNOWN_RESIDUE)


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


def backbone_trace_mmcif(mmcif_text: str) -> str:
    """A minimal mmCIF holding only the CA (protein) and P (nucleic acid) atom rows.

    A cartoon only needs the backbone trace and the viewer only needs the atom_site
    loop, so for huge entries this shrinks what is embedded in the browser iframe by
    well over an order of magnitude. It works on the text instead of re-parsing,
    which would double the cost on big files.
    """
    lines = mmcif_text.splitlines()
    try:
        first = next(i for i, line in enumerate(lines) if line.startswith("_atom_site."))
    except StopIteration:
        raise StructureParseError("mmCIF data has no atom_site loop.") from None
    rows_start = next(i for i in range(first, len(lines)) if not lines[i].startswith("_"))
    header = lines[first:rows_start]
    name_column = [line.split()[0] for line in header].index("_atom_site.label_atom_id")

    trace = ["data_trace", "loop_", *header]
    for line in lines[rows_start:]:
        if line.startswith(("#", "loop_", "_")):
            break
        if line.split()[name_column].strip("\"'") in BACKBONE_TRACE_ATOMS:
            trace.append(line)
    return "\n".join([*trace, "#"]) + "\n"


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
