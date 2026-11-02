import csv
import io

import pytest

from protein_dashboard import analysis

# Human ubiquitin, identical to the observed sequence of PDB 1UBQ chain A.
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def test_ubiquitin_is_76_residues() -> None:
    assert len(UBIQUITIN) == 76


def test_molecular_weight_matches_reference_average_mass() -> None:
    assert analysis.molecular_weight(UBIQUITIN) == pytest.approx(8564.7, abs=1.0)


def test_monoisotopic_mass_is_lighter_than_average() -> None:
    mono = analysis.molecular_weight(UBIQUITIN, monoisotopic=True)
    assert mono == pytest.approx(8559.6, abs=1.0)
    assert mono < analysis.molecular_weight(UBIQUITIN)


def test_composition_sums_to_100_and_covers_all_residues() -> None:
    composition = analysis.composition_percent(UBIQUITIN)
    assert sum(composition.values()) == pytest.approx(100.0)
    assert set(composition) == set(analysis.STANDARD_AMINO_ACIDS)
    assert composition["I"] == pytest.approx(7 / 76 * 100)


def test_property_classes_partition_the_20_amino_acids() -> None:
    grouped = "".join(analysis.PROPERTY_CLASSES.values())
    assert sorted(grouped) == sorted(analysis.STANDARD_AMINO_ACIDS)


def test_gravy_of_extreme_sequences() -> None:
    # Single-residue sequences reduce GRAVY to the raw Kyte-Doolittle value.
    assert analysis.compute_metrics("AAA").gravy == pytest.approx(1.8)
    assert analysis.compute_metrics("RRR").gravy == pytest.approx(-4.5)


def test_ubiquitin_summary_metrics() -> None:
    metrics = analysis.compute_metrics(UBIQUITIN)
    assert metrics.length == 76
    assert metrics.unknown_count == 0
    assert metrics.gravy == pytest.approx(-0.489, abs=0.001)
    # Regression values from Biopython 1.88; guard against silent scale changes.
    assert metrics.isoelectric_point == pytest.approx(6.56, abs=0.05)
    assert metrics.instability_index == pytest.approx(36.06, abs=0.05)
    assert metrics.aromaticity == pytest.approx(3 / 76, abs=1e-6)


@pytest.mark.parametrize("window", [3, 9, 19])
def test_kd_profile_length_is_sequence_minus_window_plus_one(window: int) -> None:
    profile = analysis.kd_profile(UBIQUITIN, window)
    assert len(profile.values) == len(UBIQUITIN) - window + 1
    assert len(profile.positions) == len(profile.values)


def test_kd_profile_positions_are_one_based_window_centres() -> None:
    profile = analysis.kd_profile(UBIQUITIN, 9)
    assert profile.positions[0] == 5
    assert profile.positions[-1] == 76 - 4


def test_kd_profile_flags_hydrophobic_stretch_with_tm_preset() -> None:
    sequence = "K" * 10 + "L" * 19 + "K" * 10
    profile = analysis.kd_profile(sequence, analysis.TM_WINDOW)
    above = [p for p, v in zip(profile.positions, profile.values, strict=True) if v > 1.6]
    # A window stays above 1.6 while at most 5 of its 19 residues are lysines
    # ((19-k)*3.8 - 3.9k) / 19 > 1.6  <=>  k < 5.43, so centres 15..25 qualify.
    assert above == list(range(15, 26))


@pytest.mark.parametrize("window", [1, 2, 8, 0, -3])
def test_kd_profile_rejects_even_or_tiny_windows(window: int) -> None:
    with pytest.raises(ValueError, match="odd"):
        analysis.kd_profile(UBIQUITIN, window)


def test_kd_profile_rejects_window_longer_than_sequence() -> None:
    with pytest.raises(ValueError, match="longer"):
        analysis.kd_profile("MKV", 9)


def test_contiguous_runs_end_at_breaks_and_at_non_standard_residues() -> None:
    # Break after index 2 (a gap in the density) and an unknown residue at index 6.
    assert analysis.contiguous_runs("ACDEFGXHIK", breaks=[2]) == [(0, 3), (3, 6), (7, 10)]
    assert analysis.contiguous_runs("XX") == []
    assert analysis.contiguous_runs("ACD", breaks=[2]) == [(0, 3)]  # a break after the end


def test_longest_run() -> None:
    assert analysis.longest_run("ACDEFGXHIK", breaks=[2]) == 3
    assert analysis.longest_run(UBIQUITIN) == 76
    assert analysis.longest_run("XXX") == 0


def test_kd_profile_windows_never_span_a_gap() -> None:
    # Ten leucines, a gap, ten lysines: a window across the gap would average the two.
    sequence = "L" * 10 + "K" * 10
    profile = analysis.kd_profile(sequence, 3, breaks=[9])

    assert profile.positions == [*range(2, 10), *range(12, 20)]
    # No window mixes residues across the gap: every value is pure leucine or pure lysine.
    assert profile.values == pytest.approx([3.8] * 8 + [-3.9] * 8)


def test_kd_profile_windows_never_contain_a_non_standard_residue() -> None:
    profile = analysis.kd_profile("LLLLXLLLL", 3)

    assert profile.positions == [2, 3, 7, 8]  # positions still index the full sequence
    assert profile.values == pytest.approx([3.8] * 4)


def test_kd_profile_skips_stretches_shorter_than_the_window() -> None:
    profile = analysis.kd_profile("LLLLL" + "KK", 3, breaks=[4])

    assert profile.positions == [2, 3, 4]


def test_kd_profile_rejects_window_longer_than_every_gap_free_stretch() -> None:
    with pytest.raises(ValueError, match="longer than the longest gap-free stretch \\(5"):
        analysis.kd_profile("L" * 10, 7, breaks=[4])


def test_kd_profile_of_unknown_residues_only_raises() -> None:
    with pytest.raises(ValueError, match="no standard"):
        analysis.kd_profile("XXXXX", 3)


def test_kd_profile_of_a_single_stretch_matches_biopython() -> None:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils.ProtParamData import kd

    profile = analysis.kd_profile(UBIQUITIN, 9)

    assert profile.values == pytest.approx(ProteinAnalysis(UBIQUITIN).protein_scale(kd, 9))


@pytest.mark.parametrize(
    ("length", "expected"),
    [(0, 0), (2, 0), (3, 3), (4, 3), (8, 7), (19, 19), (38, 37), (39, 39), (5000, 39)],
)
def test_max_window_is_the_largest_odd_window_that_fits(length: int, expected: int) -> None:
    assert analysis.max_window(length) == expected


def test_unknown_and_lowercase_characters_are_excluded_and_counted() -> None:
    cleaned = analysis.clean_sequence(" mk\nVX*BZ-U ")
    assert cleaned.sequence == "MKV"
    assert cleaned.unknown_count == 6  # X * B Z - U


def test_metrics_ignore_unknown_residues() -> None:
    clean = analysis.compute_metrics(UBIQUITIN)
    noisy = analysis.compute_metrics(UBIQUITIN[:10] + "XXX" + UBIQUITIN[10:])
    assert noisy.unknown_count == 3
    assert noisy.length == 76
    assert noisy.molecular_weight == pytest.approx(clean.molecular_weight)


@pytest.mark.parametrize("sequence", ["", "   ", "XXXX", "123"])
def test_sequence_without_standard_residues_raises(sequence: str) -> None:
    with pytest.raises(ValueError, match="no standard"):
        analysis.compute_metrics(sequence)
    with pytest.raises(ValueError, match="no standard"):
        analysis.molecular_weight(sequence)


def test_residue_kd_value() -> None:
    assert analysis.residue_kd_value("I") == 4.5
    assert analysis.residue_kd_value("r") == -4.5
    assert analysis.residue_kd_value("X") is None


def test_metrics_csv_roundtrip() -> None:
    metrics = analysis.compute_metrics(UBIQUITIN)
    rows = list(csv.DictReader(io.StringIO(analysis.metrics_to_csv("1UBQ_A", metrics))))
    by_name = {r["name"]: r for r in rows}
    assert float(by_name["1UBQ_A:molecular_weight_da"]["value"]) == pytest.approx(8564.74)
    assert len([r for r in rows if r["section"] == "composition_percent"]) == 20
