import pytest

from protein_dashboard import analysis, charts

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def test_composition_chart_has_one_trace_per_property_class() -> None:
    fig = charts.composition_chart(analysis.composition_percent(UBIQUITIN))

    assert [t.name for t in fig.data] == list(analysis.PROPERTY_CLASSES)
    assert sum(len(t.y) for t in fig.data) == 20
    hydrophobic = next(t for t in fig.data if t.name == "Hydrophobic")
    assert hydrophobic.y[hydrophobic.x.index("I")] == pytest.approx(7 / 76 * 100)


def test_kd_chart_hover_uses_structure_numbering_and_residue_names() -> None:
    profile = analysis.kd_profile(UBIQUITIN, 9)
    # Pretend the chain is numbered from 101, as in many PDB entries.
    numbers = list(range(101, 101 + len(UBIQUITIN)))

    fig = charts.kd_profile_chart(profile, numbers, UBIQUITIN)

    trace = fig.data[0]
    assert len(trace.y) == 76 - 9 + 1
    assert trace.x[0] == 105  # first window centre = 5th residue -> 105
    assert trace.customdata[0] == "Valine (V)"  # UBIQUITIN[4] is V
    assert "Residue %{x}" in trace.hovertemplate


def test_kd_chart_draws_threshold_line() -> None:
    profile = analysis.kd_profile(UBIQUITIN, 19)

    fig = charts.kd_profile_chart(
        profile, list(range(1, 77)), UBIQUITIN, threshold=analysis.TM_THRESHOLD
    )

    (line,) = fig.layout.shapes
    assert line.y0 == line.y1 == 1.6


def test_kd_chart_breaks_the_line_at_gaps() -> None:
    # Residues 1-10 and 21-30 are modelled; 11-20 are missing from the density.
    sequence = "L" * 10 + "K" * 10
    numbers = [*range(1, 11), *range(21, 31)]
    profile = analysis.kd_profile(sequence, 3, breaks=[9])

    trace = charts.kd_profile_chart(profile, numbers, sequence).data[0]

    gap = list(trace.y).index(None)
    assert list(trace.x[:gap]) == list(range(2, 10))
    assert trace.x[gap] is None  # Plotly leaves a hole in the line instead of joining 9 to 22
    assert list(trace.x[gap + 1 :]) == list(range(22, 30))
    assert trace.connectgaps is False


def test_kd_chart_shades_missing_residues() -> None:
    profile = analysis.kd_profile(UBIQUITIN, 9)

    fig = charts.kd_profile_chart(profile, list(range(1, 77)), UBIQUITIN, missing=[(-5, 0)])

    bands = [s for s in fig.layout.shapes if s.type == "rect"]
    assert [(b.x0, b.x1) for b in bands] == [(-5.5, 0.5)]  # covers whole residue slots


def test_kd_chart_rejects_misaligned_inputs() -> None:
    profile = analysis.kd_profile(UBIQUITIN, 9)

    with pytest.raises(ValueError, match="same length"):
        charts.kd_profile_chart(profile, [1, 2, 3], UBIQUITIN)
