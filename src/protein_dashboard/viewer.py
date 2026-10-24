"""3D view construction (py3Dmol) and the single place that renders it into Streamlit."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import py3Dmol
import streamlit.components.v1 as components

from .analysis import residue_kd_value
from .structure import ParsedStructure, trim_mmcif

# Above this many atoms the view drops surface and ligands and loads only the backbone
# trace; browsers choke on full-atom ribosome-sized entries.
HEAVY_ATOM_THRESHOLD = 50_000

KD_MIN, KD_MAX = -4.5, 4.5
HYDROPHILIC_RGB = (59, 111, 182)
NEUTRAL_RGB = (245, 245, 245)
HYDROPHOBIC_RGB = (214, 87, 29)
UNKNOWN_COLOUR = "#bdbdbd"
DIMMED_OPACITY = 0.25
LIGHT_BACKGROUND = "#ffffff"
DARK_BACKGROUND = "#0e1117"  # Streamlit's own dark page colour, so the viewer blends in


class ColorMode(StrEnum):
    CHAIN = "By chain"
    SPECTRUM = "N→C spectrum"
    HYDROPHOBICITY = "Hydrophobicity (Kyte-Doolittle)"


@dataclass(frozen=True)
class ViewerOptions:
    color_mode: ColorMode = ColorMode.CHAIN
    show_surface: bool = False
    show_ligands: bool = True
    focus_chain: str | None = None  # other chains are drawn translucent
    dark: bool = False  # the viewer lives in an iframe, so it cannot inherit the page theme


def is_heavy(structure: ParsedStructure) -> bool:
    return structure.atom_count > HEAVY_ATOM_THRESHOLD


def kd_colour(value: float) -> str:
    """Blue (hydrophilic) -> white -> orange (hydrophobic), clamped to the KD range."""
    clamped = max(KD_MIN, min(KD_MAX, value))
    if clamped < 0:
        fraction = clamped / KD_MIN
        start, end = NEUTRAL_RGB, HYDROPHILIC_RGB
    else:
        fraction = clamped / KD_MAX
        start, end = NEUTRAL_RGB, HYDROPHOBIC_RGB
    r, g, b = (round(s + (e - s) * fraction) for s, e in zip(start, end, strict=True))
    return f"#{r:02x}{g:02x}{b:02x}"


def hydrophobicity_colour_map(structure: ParsedStructure) -> dict[str, str]:
    """Residue name -> colour. Every residue of a type shares its KD value, so a per-name
    map is an exact per-residue colouring without sending one style per residue."""
    colours: dict[str, str] = {}
    for chain in structure.chains:
        for residue in chain.residues:
            value = residue_kd_value(residue.letter)
            colours[residue.name] = UNKNOWN_COLOUR if value is None else kd_colour(value)
    return colours


def _colour_style(options: ViewerOptions, structure: ParsedStructure) -> dict[str, object]:
    match options.color_mode:
        case ColorMode.CHAIN:
            return {"colorscheme": "chain"}
        case ColorMode.SPECTRUM:
            return {"color": "spectrum"}
        case ColorMode.HYDROPHOBICITY:
            return {"colorscheme": {"prop": "resn", "map": hydrophobicity_colour_map(structure)}}


def build_view(
    mmcif_text: str,
    structure: ParsedStructure,
    options: ViewerOptions,
    height: int = 600,
) -> py3Dmol.view:
    heavy = is_heavy(structure)
    view = py3Dmol.view(width="100%", height=height)
    view.setBackgroundColor(DARK_BACKGROUND if options.dark else LIGHT_BACKGROUND)
    view.addModel(trim_mmcif(mmcif_text, backbone_only=heavy), "cif")

    colours = _colour_style(options, structure)
    view.setStyle({}, {"cartoon": colours})

    others = [c.chain_id for c in structure.chains if c.chain_id != options.focus_chain]
    # With dozens of chains (ribosomes) dimming all but one would leave almost nothing visible.
    if options.focus_chain and others and not heavy:
        view.setStyle({"chain": others}, {"cartoon": {**colours, "opacity": DIMMED_OPACITY}})

    ligand_names = sorted({name for c in structure.chains for name in c.ligand_names})
    if options.show_ligands and ligand_names and not heavy:
        view.addStyle(
            {"resn": ligand_names},
            {
                "stick": {"colorscheme": "greenCarbon", "radius": 0.25},
                "sphere": {"colorscheme": "greenCarbon", "scale": 0.3},
            },
        )

    if options.show_surface and not heavy:
        surface_colour = (
            {"colorscheme": colours["colorscheme"]}
            if options.color_mode is ColorMode.HYDROPHOBICITY
            else {"color": "white"}
        )
        view.addSurface(
            py3Dmol.SAS, {"opacity": 0.6, **surface_colour}, {"resn": _polymer_names(structure)}
        )

    view.zoomTo()
    return view


def _polymer_names(structure: ParsedStructure) -> list[str]:
    return sorted({r.name for c in structure.chains for r in c.residues})


def show_view(view: py3Dmol.view, height: int = 600) -> None:
    """Render a built view in Streamlit. The only place that touches the renderer, so
    swapping py3Dmol's HTML export for another component later is a one-function change."""
    components.html(view._make_html(), height=height)
