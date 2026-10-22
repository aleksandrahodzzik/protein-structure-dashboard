from pathlib import Path

import pytest

from protein_dashboard import structure
from protein_dashboard.structure import StructureParseError, parse_structure

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
DATA_DIR = Path(__file__).parent / "data"

ATOM_SITE_HEADER = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""

BACKBONE = ("N", "CA", "C")


def build_mmcif(
    records: list[tuple[str, str, str, tuple[str, ...]]],
    gaps_after: frozenset[int] = frozenset(),
) -> str:
    """records: (group, residue name, chain, atom names). Residues are laid out along x,
    3.8 A apart, so consecutive backbones are 1.3 A C->N apart and PPBuilder links them.
    Residues after an index in ``gaps_after`` are shifted 10 A away, i.e. a chain break."""
    lines = [ATOM_SITE_HEADER.rstrip("\n")]
    serial = 1
    for index, (group, resname, chain, atoms) in enumerate(records, start=1):
        for offset, atom in enumerate(atoms):
            shift = 10.0 * sum(gap < index for gap in gaps_after)
            x = 3.8 * index + shift + 1.2 * offset
            y = 0.0 if atom in BACKBONE else 2.0 + offset
            element = atom[0] if atom != "CA" or resname != "CA" else "CA"
            lines.append(
                f"{group} {serial} {element} {atom} . {resname} {chain} 1 {index} ? "
                f"{x:.3f} {y:.3f} 0.000 1.00 10.00 ? {index} {resname} {chain} {atom} 1"
            )
            serial += 1
    return "\n".join(lines) + "\n#\n"


def residue(
    resname: str, chain: str = "A", group: str = "ATOM"
) -> tuple[str, str, str, tuple[str, ...]]:
    return (group, resname, chain, BACKBONE)


@pytest.fixture(scope="module")
def ubiquitin() -> structure.ParsedStructure:
    text = (DATA_DIR / "1UBQ.cif").read_text()
    return parse_structure(text, pdb_id="1UBQ")


def test_1ubq_has_one_protein_chain_of_76_residues(ubiquitin: structure.ParsedStructure) -> None:
    assert [c.chain_id for c in ubiquitin.chains] == ["A"]
    chain = ubiquitin.chain("A")
    assert chain.residue_count == 76
    assert chain.sequence == UBIQUITIN
    assert chain.unknown_count == 0
    assert chain.segment_count == 1


def test_1ubq_waters_are_not_ligands_and_are_counted_as_atoms(
    ubiquitin: structure.ParsedStructure,
) -> None:
    assert ubiquitin.chain("A").ligands == ()
    assert ubiquitin.atom_count == 660  # 602 protein atoms + 58 waters


def test_1ubq_residue_records_carry_numbering(ubiquitin: structure.ParsedStructure) -> None:
    first, last = ubiquitin.chain("A").residues[0], ubiquitin.chain("A").residues[-1]
    assert (first.number, first.name, first.letter) == (1, "MET", "M")
    assert (last.number, last.name, last.letter) == (76, "GLY", "G")


def test_unknown_chain_id_raises_key_error(ubiquitin: structure.ParsedStructure) -> None:
    with pytest.raises(KeyError):
        ubiquitin.chain("Z")


@pytest.mark.parametrize(
    ("resname", "letter"),
    [
        ("ALA", "A"),
        ("ala", "A"),
        ("MSE", "M"),  # selenomethionine
        ("SEP", "S"),  # phosphoserine
        ("TPO", "T"),  # phosphothreonine
        ("PTR", "Y"),  # phosphotyrosine
        ("UNK", "X"),
        ("ZZZ", "X"),
        ("HOH", "X"),
    ],
)
def test_residue_letter_mapping(resname: str, letter: str) -> None:
    assert structure.residue_letter(resname) == letter


def test_modified_and_unknown_residues_stay_in_sequence_in_order() -> None:
    mmcif = build_mmcif(
        [
            residue("ALA"),
            residue("MSE", group="HETATM"),
            residue("UNK"),
            residue("GLY"),
        ]
    )

    chain = parse_structure(mmcif, "TEST").chain("A")

    assert chain.sequence == "AMXG"
    assert chain.residue_count == 4
    assert chain.unknown_count == 1
    assert [r.name for r in chain.residues] == ["ALA", "MSE", "UNK", "GLY"]
    # Modified residues written as HETATM are polymer residues, not ligands.
    assert chain.ligands == ()


def test_isolated_residues_between_gaps_are_kept_and_counted_as_segments() -> None:
    # A G | K | S T: the lone lysine has no bonded neighbour, which PPBuilder alone would drop.
    mmcif = build_mmcif(
        [residue("ALA"), residue("GLY"), residue("LYS"), residue("SER"), residue("THR")],
        gaps_after=frozenset({2, 3}),
    )

    chain = parse_structure(mmcif, "TEST").chain("A")

    assert chain.sequence == "AGKST"
    assert chain.segment_count == 3


def test_ligands_exclude_water_and_ions_do_not_pollute_the_sequence() -> None:
    mmcif = build_mmcif(
        [
            residue("ALA"),
            residue("GLY"),
            ("HETATM", "CA", "A", ("CA",)),  # calcium ion: residue and atom both named CA
            ("HETATM", "HEM", "A", ("FE", "NA", "NB")),
            ("HETATM", "HOH", "A", ("O",)),
        ]
    )

    chain = parse_structure(mmcif, "TEST").chain("A")

    assert chain.sequence == "AG"
    assert chain.ligand_names == ("CA", "HEM")
    assert len(chain.ligands) == 2


def test_multiple_chains_and_nucleic_acid_chains_are_separated() -> None:
    mmcif = build_mmcif(
        [
            residue("ALA", "A"),
            residue("GLY", "A"),
            residue("LYS", "B"),
            ("ATOM", "DA", "C", ("P", "OP1", "O5'")),  # DNA-only chain
        ]
    )

    parsed = parse_structure(mmcif, "TEST")

    assert [c.chain_id for c in parsed.chains] == ["A", "B"]
    assert parsed.chain("B").sequence == "K"


def test_structure_without_protein_is_rejected() -> None:
    mmcif = build_mmcif([("ATOM", "DA", "A", ("P", "OP1", "O5'"))])

    with pytest.raises(StructureParseError, match="no protein"):
        parse_structure(mmcif, "TEST")


@pytest.mark.parametrize("text", ["", "not a cif file", "data_X\n#\n"])
def test_garbage_input_raises_parse_error(text: str) -> None:
    with pytest.raises(StructureParseError):
        parse_structure(text)
