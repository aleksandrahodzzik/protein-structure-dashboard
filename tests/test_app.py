"""Smoke tests: drive the real Streamlit script with a fake RCSB client (no network)."""

from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path

import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

from helpers import build_mmcif, poly_seq_scheme, residue
from protein_dashboard import app, examples
from protein_dashboard.rcsb_client import (
    EntryMetadata,
    EntryNotFoundError,
    RcsbUnavailableError,
)

ROOT = Path(__file__).parents[1]
UBIQUITIN_CIF = (ROOT / "tests" / "data" / "1UBQ.cif").read_text()
UBIQUITIN = EntryMetadata(
    pdb_id="1UBQ",
    title="STRUCTURE OF UBIQUITIN REFINED AT 1.8 ANGSTROMS RESOLUTION",
    method="X-RAY DIFFRACTION",
    resolution=1.8,
    organisms=("Homo sapiens",),
    deposition_date="1987-01-02",
    atom_count=660,
)


# A synthetic chain with every way the observed sequence can differ from SEQRES:
# residues 13-20 are missing (a gap), residue 32 is LYS where SEQRES says ARG, and
# residue 99 is modelled but has no SEQRES position.
GAPPED_NUMBERS = [*range(1, 13), *range(21, 33), 99]
GAPPED_CIF = build_mmcif(
    [residue(name) for name in ["LEU"] * 12 + ["LYS"] * 12 + ["SER"]],
    gaps_after=frozenset({12}),
    numbering=GAPPED_NUMBERS,
) + poly_seq_scheme(
    [("A", n, "LEU") for n in range(1, 13)]
    + [("A", n, "GLY") for n in range(13, 21)]
    + [("A", n, "LYS") for n in range(21, 32)]
    + [("A", 32, "ARG")]
)
NO_SEQRES_CIF = build_mmcif([residue("LEU")] * 10)


def synthetic_entry(pdb_id: str) -> EntryMetadata:
    return replace(UBIQUITIN, pdb_id=pdb_id, title=f"SYNTHETIC {pdb_id}", atom_count=None)


ENTRIES = {
    "1UBQ": (UBIQUITIN, UBIQUITIN_CIF),
    "2GAP": (synthetic_entry("2GAP"), GAPPED_CIF),
    "3NOS": (synthetic_entry("3NOS"), NO_SEQRES_CIF),
}


class FakeClient:
    """Stands in for RcsbClient; only ENTRIES exist, everything else is unknown."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch_entry(self, pdb_id: str, *, with_organisms: bool = True) -> EntryMetadata:
        self.calls.append(f"entry:{pdb_id}")
        if pdb_id not in ENTRIES:
            raise EntryNotFoundError(f"RCSB has no entry for '{pdb_id}'.")
        return ENTRIES[pdb_id][0]

    def fetch_structure(self, pdb_id: str) -> str:
        self.calls.append(f"structure:{pdb_id}")
        if pdb_id not in ENTRIES:
            raise EntryNotFoundError(f"RCSB has no entry for '{pdb_id}'.")
        return ENTRIES[pdb_id][1]

    def search(self, query: str, limit: int = 10) -> list[str]:
        self.calls.append(f"search:{query}")
        return ["1UBQ"] if "ubiquitin" in query.lower() else []


@pytest.fixture
def client(monkeypatch: pytest.MonkeyPatch) -> Iterator[FakeClient]:
    fake = FakeClient()
    monkeypatch.setattr(app, "get_client", lambda: fake)
    st.cache_data.clear()  # cached results from another test must not leak in
    yield fake
    st.cache_data.clear()


def new_app() -> AppTest:
    return AppTest.from_file(str(ROOT / "streamlit_app.py"), default_timeout=60)


def submit_search(at: AppTest, query: str) -> AppTest:
    at.sidebar.text_input[0].set_value(query)
    next(b for b in at.sidebar.button if b.label == "Search").click()
    return at.run()


def metric_values(at: AppTest) -> dict[str, str]:
    return {m.label: m.value for m in at.metric}


def card_labels(at: AppTest) -> list[str]:
    return [b.label for b in at.main.button]


def test_landing_page_starts_on_the_first_category(client: FakeClient) -> None:
    at = new_app().run()

    assert not at.exception
    assert at.title[0].value == "Protein Structure Dashboard"
    assert card_labels(at) == [e.pdb_id for e in examples.CATEGORIES[0].items]
    assert at.pills[0].options == [c.name for c in examples.CATEGORIES]
    assert client.calls == []  # nothing is fetched until the user chooses


def test_every_category_shows_exactly_its_own_examples(client: FakeClient) -> None:
    at = new_app().run()

    for category in examples.CATEGORIES:
        at.pills[0].set_value(category.name)
        at.run()

        assert not at.exception, category.name
        assert card_labels(at) == [e.pdb_id for e in category.items], category.name


def test_home_button_returns_from_an_entry_to_the_example_gallery(client: FakeClient) -> None:
    at = new_app()
    at.query_params["pdb"] = "1UBQ"
    at.run()
    assert metric_values(at)  # an entry is open

    at.sidebar.button(key="home").click()
    at.run()

    assert not at.exception
    assert "pdb" not in at.query_params
    assert card_labels(at) == [e.pdb_id for e in examples.CATEGORIES[0].items]


def theme_scripts(at: AppTest) -> list[str]:
    return [h.value for h in at.html if "stMainMenuItem-theme" in h.value]


def test_theme_button_requests_the_opposite_theme_once(client: FakeClient) -> None:
    at = new_app().run()
    assert theme_scripts(at) == []

    at.sidebar.button(key="theme_toggle").click()
    at.run()

    assert not at.exception
    assert "theme_target" not in at.session_state  # consumed by the script that applies it
    (script,) = theme_scripts(at)
    assert "stMainMenuItem-theme-Dark" in script  # AppTest has no theme, i.e. light

    at.run()
    assert theme_scripts(at) == []  # sent once, not on every rerun


def theme_icon(at: AppTest) -> str:
    return at.sidebar.button(key="theme_toggle").label


def test_theme_button_shows_the_new_theme_without_a_rerun(client: FakeClient) -> None:
    # The menu that viewers of a deployed app get has no "Rerun", so the icon must not wait
    # for st.context.theme to catch up (AppTest never reports a theme, i.e. light).
    at = new_app().run()
    assert theme_icon(at) == ":material/dark_mode:"

    at.sidebar.button(key="theme_toggle").click()
    at.run()

    assert theme_icon(at) == ":material/light_mode:"
    (script,) = theme_scripts(at)
    assert "rerun" not in script

    at.sidebar.button(key="theme_toggle").click()
    at.run()

    assert theme_icon(at) == ":material/dark_mode:"
    (script,) = theme_scripts(at)
    assert "stMainMenuItem-theme-Light" in script


class FakeContext:
    def __init__(self, theme: str | None) -> None:
        self.theme = {"type": theme}


def test_theme_follows_the_browser_once_it_reports_a_change(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    state: dict[str, bool] = {}
    monkeypatch.setattr(app.st, "session_state", state)
    monkeypatch.setattr(app.st, "context", FakeContext("light"))
    assert not app.is_dark_theme()

    app.request_theme_switch()  # the button asked for dark, the browser still says light
    assert app.is_dark_theme()
    assert state["theme_target"] == "Dark"

    monkeypatch.setattr(app.st, "context", FakeContext("dark"))  # the switch went through
    assert app.is_dark_theme()

    monkeypatch.setattr(app.st, "context", FakeContext("light"))  # Streamlit's own menu
    assert not app.is_dark_theme()


def open_entry(pdb_id: str) -> AppTest:
    at = new_app()
    at.query_params["pdb"] = pdb_id
    return at.run()


def captions(at: AppTest) -> str:
    return "\n".join(c.value for c in at.caption)


def test_observed_sequence_matches_seqres_for_1ubq(client: FakeClient) -> None:
    at = open_entry("1UBQ")

    assert not at.exception
    assert "76 of 76 SEQRES residues are modelled (100%), in 1 contiguous segment." in captions(at)
    assert "Not modelled" not in captions(at)
    assert "line breaks" not in captions(at)  # nothing to explain under the profile


def test_differences_from_seqres_are_shown(client: FakeClient) -> None:
    at = open_entry("2GAP")

    assert not at.exception
    text = captions(at)
    assert "24 of 32 SEQRES residues are modelled (75%), in 2 contiguous segments" in text
    assert "Not modelled (8 residues): 13\u201320." in text
    assert at.warning[0].value == (
        "1 modelled residue(s) differ from SEQRES: 32 (ARG in SEQRES, LYS modelled)"
    )
    assert "have no position in SEQRES and are still analysed: 99." in at.info[0].value
    assert at.expander[0].label == "Observed sequence against SEQRES"
    assert at.code[0].value.splitlines()[:2] == [
        "SEQRES 1 LLLLLLLLLL LLGGGGGGGG KKKKKKKKKK KR",
        "Model    LLLLLLLLLL LL-------- KKKKKKKKKK KK",
    ]


def test_profile_breaks_at_the_gap_and_the_window_fits_the_longest_stretch(
    client: FakeClient,
) -> None:
    at = open_entry("2GAP")

    # 12 residues before the gap, 13 after it (21-32 and 99): the largest odd window is 13.
    assert at.sidebar.slider[0].max == 13
    assert "the line breaks there" in captions(at)
    assert "Grey bands mark SEQRES residues that are not modelled." in captions(at)


def test_file_without_seqres_says_so(client: FakeClient) -> None:
    at = open_entry("3NOS")

    assert not at.exception
    assert "This file has no SEQRES mapping for the chain" in captions(at)
    assert not at.expander


def test_entering_1ubq_renders_metadata_and_metrics(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "1ubq")

    assert not at.exception
    assert at.query_params["pdb"] == "1UBQ"  # the selection lives in the URL, so it is shareable
    assert "STRUCTURE OF UBIQUITIN" in at.subheader[0].value
    metrics = metric_values(at)
    assert metrics["Length"] == "76 aa"
    assert metrics["Molecular weight"] == "8,564.7 Da"
    assert metrics["Isoelectric point (pI)"] == "6.56"
    assert metrics["GRAVY"] == "-0.489"


def test_clicking_an_example_loads_it(client: FakeClient) -> None:
    at = new_app().run()

    at.main.button[0].click()
    at.run()

    assert not at.exception
    assert metric_values(at)["Length"] == "76 aa"


def test_second_view_of_the_same_entry_is_served_from_cache(client: FakeClient) -> None:
    submit_search(new_app().run(), "1UBQ")
    network_calls = len(client.calls)

    submit_search(new_app().run(), "1UBQ")

    assert len(client.calls) == network_calls


def test_transmembrane_preset_sets_window_19(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "1UBQ")
    assert at.sidebar.slider[0].value == 9

    preset = next(b for b in at.sidebar.button if b.label.startswith("Transmembrane preset"))
    preset.click()
    at.run()

    assert not at.exception
    assert at.sidebar.slider[0].value == 19


def test_slider_range_is_valid_odd_windows(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "1UBQ")

    slider = at.sidebar.slider[0]
    assert (slider.min, slider.max, slider.step) == (3, 39, 2)


def test_unknown_id_shows_a_clear_error(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "0ZZZ")

    assert not at.exception
    assert at.error[0].value == "RCSB has no entry for '0ZZZ'."
    assert not at.metric


def test_rcsb_outage_is_reported_not_raised(
    client: FakeClient, monkeypatch: pytest.MonkeyPatch
) -> None:
    def down(pdb_id: str, **_: object) -> EntryMetadata:
        raise RcsbUnavailableError("RCSB did not respond successfully after 4 attempts.")

    monkeypatch.setattr(client, "fetch_entry", down)

    at = submit_search(new_app().run(), "1UBQ")

    assert not at.exception
    assert "did not respond" in at.error[0].value


def test_keyword_search_lists_results_and_opens_one(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "ubiquitin")

    assert not at.exception
    assert client.calls[0] == "search:ubiquitin"
    result = next(b for b in at.sidebar.button if b.label.startswith("1UBQ"))
    result.click()
    at.run()

    assert metric_values(at)["Length"] == "76 aa"


def test_keyword_search_without_hits_says_so(client: FakeClient) -> None:
    at = submit_search(new_app().run(), "nonsense")

    assert not at.exception
    assert "No structures found for 'nonsense'" in at.sidebar.warning[0].value
