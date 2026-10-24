from dataclasses import replace
from pathlib import Path

import pytest

from helpers import build_mmcif, residue
from protein_dashboard import viewer
from protein_dashboard.structure import (
    ParsedStructure,
    StructureParseError,
    parse_structure,
    trim_mmcif,
)
from protein_dashboard.viewer import ColorMode, ViewerOptions, build_view

CIF = (Path(__file__).parent / "data" / "1UBQ.cif").read_text()


@pytest.fixture(scope="module")
def ubiquitin() -> ParsedStructure:
    return parse_structure(CIF, "1UBQ")


def html_for(structure: ParsedStructure, options: ViewerOptions, text: str = CIF) -> str:
    return build_view(text, structure, options)._make_html()


def test_kd_colour_endpoints_and_clamping() -> None:
    assert viewer.kd_colour(0) == "#f5f5f5"
    assert viewer.kd_colour(4.5) == "#d6571d"
    assert viewer.kd_colour(-4.5) == "#3b6fb6"
    assert viewer.kd_colour(100) == viewer.kd_colour(4.5)
    assert viewer.kd_colour(-100) == viewer.kd_colour(-4.5)


def test_hydrophobicity_map_links_residue_names_to_kd_values(ubiquitin: ParsedStructure) -> None:
    colours = viewer.hydrophobicity_colour_map(ubiquitin)

    assert colours["ILE"] == viewer.kd_colour(4.5)
    assert colours["ARG"] == viewer.kd_colour(-4.5)
    assert set(colours) == {r.name for r in ubiquitin.chain("A").residues}


def test_modified_residue_gets_its_parent_residue_colour() -> None:
    mmcif = build_mmcif([residue("ALA"), residue("MSE", group="HETATM"), residue("UNK")])
    colours = viewer.hydrophobicity_colour_map(parse_structure(mmcif, "T"))

    assert colours["MSE"] == viewer.kd_colour(1.9)  # methionine
    assert colours["UNK"] == viewer.UNKNOWN_COLOUR


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        (ColorMode.CHAIN, '"colorscheme":"chain"'),
        (ColorMode.SPECTRUM, '"color":"spectrum"'),
        (ColorMode.HYDROPHOBICITY, '"prop":"resn"'),
    ],
)
def test_each_colour_mode_emits_its_style(
    ubiquitin: ParsedStructure, mode: ColorMode, expected: str
) -> None:
    html = html_for(ubiquitin, ViewerOptions(color_mode=mode))

    assert '"cartoon"' in html
    assert expected in html.replace(" ", "")


def test_surface_is_off_by_default_and_toggles_on(ubiquitin: ParsedStructure) -> None:
    assert "addSurface" not in html_for(ubiquitin, ViewerOptions())
    assert "addSurface" in html_for(ubiquitin, ViewerOptions(show_surface=True))


def test_ligands_are_drawn_as_sticks_when_present() -> None:
    mmcif = build_mmcif(
        [residue("ALA"), residue("GLY"), ("HETATM", "HEM", "A", ("FE", "NA", "NB"))]
    )
    parsed = parse_structure(mmcif, "T")

    on = html_for(parsed, ViewerOptions(show_ligands=True), mmcif)
    off = html_for(parsed, ViewerOptions(show_ligands=False), mmcif)

    assert '"stick"' in on
    assert '"HEM"' in on
    assert '"stick"' not in off


def test_focus_chain_dims_the_other_chains() -> None:
    mmcif = build_mmcif([residue("ALA", "A"), residue("GLY", "A"), residue("LYS", "B")])
    parsed = parse_structure(mmcif, "T")

    html = html_for(parsed, ViewerOptions(focus_chain="A"), mmcif)

    assert '"chain":["B"]' in html.replace(" ", "")
    assert '"opacity":0.25' in html.replace(" ", "")


def test_heavy_multi_chain_structures_are_not_dimmed() -> None:
    mmcif = build_mmcif([residue("ALA", "A"), residue("GLY", "A"), residue("LYS", "B")])
    heavy = replace(parse_structure(mmcif, "T"), atom_count=viewer.HEAVY_ATOM_THRESHOLD + 1)

    html = html_for(heavy, ViewerOptions(focus_chain="A"), mmcif)

    assert '"opacity":0.25' not in html.replace(" ", "")


def test_heavy_structures_get_backbone_only_and_no_surface_or_ligands_or_dimming(
    ubiquitin: ParsedStructure,
) -> None:
    heavy = replace(ubiquitin, atom_count=viewer.HEAVY_ATOM_THRESHOLD + 1)
    options = ViewerOptions(show_surface=True, show_ligands=True, focus_chain="A")

    html = html_for(heavy, options)

    assert viewer.is_heavy(heavy)
    assert not viewer.is_heavy(ubiquitin)
    assert "addSurface" not in html
    # Atom rows are embedded as one JS string, so rows are separated by a literal "\n".
    assert html.count("\\nATOM ") == 76 * 4  # backbone atoms only
    assert html_for(ubiquitin, options).count("\\nATOM ") == 602


def test_backbone_only_trim_keeps_cartoon_atoms_and_stays_parseable() -> None:
    trimmed = trim_mmcif(CIF, backbone_only=True)

    parsed = parse_structure(trimmed, "1UBQ")
    assert parsed.atom_count == 76 * 4  # N, CA, C, O per residue; side chains and waters gone
    assert len(trimmed) < len(CIF) / 2
    assert "_struct_conf." in trimmed  # helix/sheet records survive so the cartoon keeps its shape
    assert "_citation." not in trimmed


def test_trim_without_backbone_only_keeps_every_atom_of_the_first_model() -> None:
    assert parse_structure(trim_mmcif(CIF), "1UBQ").atom_count == 660


def test_trim_drops_later_models_of_an_nmr_ensemble() -> None:
    ensemble = build_mmcif([residue("ALA"), residue("GLY")], models=3)
    assert "_atom_site.pdbx_PDB_model_num" in ensemble

    trimmed = trim_mmcif(ensemble)

    assert trimmed.count("ATOM ") == 6  # 2 residues x 3 backbone atoms, model 1 only
    assert " 2\n" not in trimmed.replace(" 2 ", "")  # no rows tagged with model 2 or 3
    assert parse_structure(trimmed, "T").atom_count == 6


def test_trim_rejects_text_without_atoms() -> None:
    with pytest.raises(StructureParseError, match="atom_site"):
        trim_mmcif("data_X\n#\n")


def test_show_view_renders_through_components_html(
    ubiquitin: ParsedStructure, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[str, int]] = []
    monkeypatch.setattr(
        viewer.components, "html", lambda html, height: calls.append((html, height))
    )

    viewer.show_view(build_view(CIF, ubiquitin, ViewerOptions()), height=480)

    assert calls[0][1] == 480
    assert "3dmolviewer" in calls[0][0]
