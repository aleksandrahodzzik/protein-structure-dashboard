"""Plotly figures. They take plain data from `analysis` and know nothing about Streamlit."""

from __future__ import annotations

from collections.abc import Sequence

import plotly.graph_objects as go

from .analysis import (
    AMINO_ACID_NAMES,
    PROPERTY_CLASSES,
    TM_THRESHOLD,
    HydrophobicityProfile,
)

CLASS_COLOURS = {
    "Hydrophobic": "#d95f02",
    "Polar": "#1b9e77",
    "Charged": "#7570b3",
    "Special": "#666666",
}
PROFILE_COLOUR = "#1f77b4"


def composition_chart(composition: dict[str, float]) -> go.Figure:
    """Residue percentages, one trace per property class so the legend explains colour."""
    fig = go.Figure()
    for group, letters in PROPERTY_CLASSES.items():
        fig.add_bar(
            x=list(letters),
            y=[composition[aa] for aa in letters],
            name=group,
            marker_color=CLASS_COLOURS[group],
            customdata=[AMINO_ACID_NAMES[aa] for aa in letters],
            hovertemplate="%{customdata} (%{x}): %{y:.1f}%<extra>" + group + "</extra>",
        )
    fig.update_layout(
        title="Amino acid composition",
        xaxis={
            "title": "Residue, grouped by property class",
            "categoryorder": "array",
            "categoryarray": [aa for letters in PROPERTY_CLASSES.values() for aa in letters],
        },
        yaxis={"title": "Share of residues (%)"},
        legend={"orientation": "h", "y": -0.25},
        margin={"t": 50, "b": 40},
    )
    return fig


def kd_profile_chart(
    profile: HydrophobicityProfile,
    residue_numbers: Sequence[int],
    sequence: str,
    threshold: float = TM_THRESHOLD,
) -> go.Figure:
    """Kyte-Doolittle profile; hover names the window-centre residue.

    ``sequence`` and ``residue_numbers`` describe the residues the profile was computed
    from (same length, same order), so hover shows the structure's own numbering.
    """
    if len(residue_numbers) != len(sequence):
        raise ValueError("residue_numbers and sequence must have the same length.")
    centres = [position - 1 for position in profile.positions]  # 0-based indices
    x = [residue_numbers[c] for c in centres]
    names = [f"{AMINO_ACID_NAMES[sequence[c]]} ({sequence[c]})" for c in centres]

    fig = go.Figure()
    fig.add_scatter(
        x=x,
        y=profile.values,
        mode="lines",
        line={"color": PROFILE_COLOUR, "width": 2},
        name=f"Window {profile.window}",
        customdata=names,
        hovertemplate="Residue %{x} %{customdata}<br>KD score: %{y:.2f}<extra></extra>",
    )
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="#b2182b",
        annotation_text=f"Threshold {threshold:g} (TM helix, window 19)",
        annotation_position="top left",
    )
    fig.update_layout(
        title=f"Kyte-Doolittle hydrophobicity (window {profile.window})",
        xaxis={"title": "Residue number (window centre)"},
        yaxis={"title": "Mean hydropathy"},
        showlegend=False,
        margin={"t": 50, "b": 40},
    )
    return fig
