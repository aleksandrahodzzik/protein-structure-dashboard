from pathlib import Path

import pytest

from helpers import build_mmcif, poly_seq_scheme, residue
from protein_dashboard import structure
from protein_dashboard.structure import StructureParseError, parse_structure

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
DATA_DIR = Path(__file__).parent / "data"


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


def test_residue_letter_always_returns_a_standard_letter_or_x() -> None:
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended

    allowed = set(structure.STANDARD_AMINO_ACIDS) | {structure.UNKNOWN_RESIDUE}
    assert {structure.residue_letter(name) for name in protein_letters_3to1_extended} <= allowed


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
    assert [r.segment for r in chain.residues] == [0, 0, 1, 2, 2]
    assert chain.breaks == (1, 2)  # gaps after the 2nd and the 3rd residue


def test_unbroken_chain_has_no_breaks(ubiquitin: structure.ParsedStructure) -> None:
    assert ubiquitin.chain("A").breaks == ()


def test_ca_only_models_are_linked_by_ca_distance() -> None:
    # Some low-resolution models have only CA atoms: there is no C-N bond to measure,
    # so consecutive residues 3.8 A apart must still count as neighbours.
    ca_only = [("ATOM", name, "A", ("CA",)) for name in ("ALA", "GLY", "LYS", "SER")]

    chain = parse_structure(build_mmcif(ca_only, gaps_after=frozenset({2})), "T").chain("A")

    assert chain.sequence == "AGKS"
    assert chain.breaks == (1,)


def test_residue_without_backbone_atoms_cannot_be_linked() -> None:
    # A standard residue modelled with its side chain only: no C, N or CA to measure.
    records = [residue("ALA"), ("ATOM", "LEU", "A", ("CB", "CG")), residue("GLY")]

    chain = parse_structure(build_mmcif(records), "T").chain("A")

    assert chain.sequence == "ALG"
    assert chain.breaks == (0, 1)


def test_seqres_rows_without_an_author_number_are_skipped() -> None:
    text = build_mmcif([residue("ALA")]) + poly_seq_scheme([("A", 1, "ALA"), ("A", 2, "GLY")])

    chain = parse_structure(text.replace(" 2 ? ? ? A . n", " ? ? ? ? A . n"), "T").chain("A")

    assert [s.number for s in chain.seqres] == [1]


def test_malformed_seqres_block_is_ignored() -> None:
    text = build_mmcif([residue("ALA")]) + poly_seq_scheme([("A", 1, "ALA")])

    chain = parse_structure(text.replace("_pdbx_poly_seq_scheme.pdb_seq_num\n", ""), "T")

    assert chain.chain("A").seqres == ()


def test_seqres_loop_with_a_missing_value_does_not_break_parsing() -> None:
    text = build_mmcif([residue("ALA")]) + poly_seq_scheme([("A", 1, "ALA"), ("A", 2, "GLY")])
    # Make pdb_ins_code the last column, then drop its value from the last row: the
    # column ends up one value short of the others.
    text = text.replace("_pdbx_poly_seq_scheme.hetero\n", "").replace(" . n\n", " .\n")
    text = text.replace("A .\n#", "A\n#")

    parsed = parse_structure(text, "T")

    assert parsed.chain("A").sequence == "A"
    assert parsed.chain("A").seqres == ()


def test_seqres_is_read_per_author_chain(ubiquitin: structure.ParsedStructure) -> None:
    seqres = ubiquitin.chain("A").seqres

    assert "".join(s.letter for s in seqres) == UBIQUITIN
    assert (seqres[0].number, seqres[0].name) == (1, "MET")


def test_1ubq_observed_sequence_matches_seqres(ubiquitin: structure.ParsedStructure) -> None:
    comparison = structure.compare_to_seqres(ubiquitin.chain("A"))

    assert comparison is not None
    assert (comparison.seqres_length, comparison.modelled_count) == (76, 76)
    assert comparison.missing == ()
    assert comparison.mismatches == ()
    assert comparison.unmatched == ()


def gapped_chain(seqres_rows: list[tuple[str, int, str]]) -> structure.ProteinChain:
    """Modelled: ALA 2, GLY 3 | gap | LYS 6. SEQRES is supplied by the test."""
    mmcif = build_mmcif(
        [residue("ALA"), residue("GLY"), residue("LYS")],
        gaps_after=frozenset({2}),
        numbering=[2, 3, 6],
    )
    return parse_structure(mmcif + poly_seq_scheme(seqres_rows), "T").chain("A")


SEQRES_MAGSTK = [
    ("A", n, name) for n, name in enumerate(["MET", "ALA", "GLY", "SER", "THR", "LYS"], 1)
]


def test_missing_residues_are_listed_as_stretches_in_author_numbering() -> None:
    comparison = structure.compare_to_seqres(gapped_chain(SEQRES_MAGSTK))

    assert comparison is not None
    assert (comparison.seqres_length, comparison.modelled_count) == (6, 3)
    assert [(m.label, m.count) for m in comparison.missing] == [("1", 1), ("4\u20135", 2)]
    assert comparison.missing_count == 3
    assert comparison.mismatches == ()


def test_missing_stretch_from_a_negative_number_is_spelled_out() -> None:
    # Expression tags are often numbered below 1 (e.g. -18 to 26 in 6VXX).
    rows = [("A", n, "GLY") for n in (-1, 0, 1)] + [("A", 2, "ALA"), ("A", 3, "GLY")]

    comparison = structure.compare_to_seqres(gapped_chain(rows))

    assert comparison is not None
    assert [m.label for m in comparison.missing] == ["-1 to 1"]


def test_alignment_marks_missing_residues_under_the_seqres_line() -> None:
    comparison = structure.compare_to_seqres(gapped_chain(SEQRES_MAGSTK))
    assert comparison is not None

    assert comparison.alignment(width=4, block=2) == (
        "SEQRES 1 MA GS\nModel    -A G-\n\nSEQRES 5 TK\nModel    -K"
    )


def test_modelled_residue_that_differs_from_seqres_is_a_mismatch() -> None:
    rows = [("A", 2, "ALA"), ("A", 3, "SER"), ("A", 6, "LYS")]  # SEQRES says SER at 3

    comparison = structure.compare_to_seqres(gapped_chain(rows))

    assert comparison is not None
    assert [(m.label, m.seqres_name, m.modelled_name) for m in comparison.mismatches] == [
        ("3", "SER", "GLY")
    ]


def test_heterogeneous_seqres_position_accepts_either_residue() -> None:
    rows = [("A", 2, "ALA"), ("A", 3, "SER"), ("A", 3, "GLY"), ("A", 6, "LYS")]

    comparison = structure.compare_to_seqres(gapped_chain(rows))

    assert comparison is not None
    assert comparison.seqres_length == 3
    assert comparison.mismatches == ()


def test_modelled_residue_outside_seqres_is_reported_as_unmatched() -> None:
    rows = [("A", 2, "ALA"), ("A", 3, "GLY")]  # residue 6 is not in the deposited sequence

    comparison = structure.compare_to_seqres(gapped_chain(rows))

    assert comparison is not None
    assert [r.number for r in comparison.unmatched] == [6]


def test_seqres_of_other_chains_is_not_mixed_in() -> None:
    rows = [("B", n, "GLY") for n in range(1, 10)]

    assert structure.compare_to_seqres(gapped_chain(rows)) is None


def test_file_without_seqres_mapping_gives_no_comparison() -> None:
    chain = parse_structure(build_mmcif([residue("ALA"), residue("GLY")]), "T").chain("A")

    assert chain.seqres == ()
    assert structure.compare_to_seqres(chain) is None


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
