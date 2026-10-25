"""Smoke tests: drive the real Streamlit script with a fake RCSB client (no network)."""

from collections.abc import Iterator
from pathlib import Path

import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

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


class FakeClient:
    """Stands in for RcsbClient; only 1UBQ exists, everything else is unknown."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch_entry(self, pdb_id: str, *, with_organisms: bool = True) -> EntryMetadata:
        self.calls.append(f"entry:{pdb_id}")
        if pdb_id != "1UBQ":
            raise EntryNotFoundError(f"RCSB has no entry for '{pdb_id}'.")
        return UBIQUITIN

    def fetch_structure(self, pdb_id: str) -> str:
        self.calls.append(f"structure:{pdb_id}")
        if pdb_id != "1UBQ":
            raise EntryNotFoundError(f"RCSB has no entry for '{pdb_id}'.")
        return UBIQUITIN_CIF

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


def test_landing_page_offers_every_example_entry(client: FakeClient) -> None:
    at = new_app().run()

    assert not at.exception
    assert at.title[0].value == "Protein Structure Dashboard"
    labels = [b.label for b in at.main.button]
    assert labels == [e.pdb_id for e in examples.ALL]
    assert client.calls == []  # nothing is fetched until the user chooses


def test_home_button_returns_from_an_entry_to_the_example_gallery(client: FakeClient) -> None:
    at = new_app()
    at.query_params["pdb"] = "1UBQ"
    at.run()
    assert metric_values(at)  # an entry is open

    at.sidebar.button(key="home").click()
    at.run()

    assert not at.exception
    assert "pdb" not in at.query_params
    assert [b.label for b in at.main.button] == [e.pdb_id for e in examples.ALL]


def test_theme_button_requests_the_opposite_theme_once(client: FakeClient) -> None:
    at = new_app().run()

    at.sidebar.button(key="theme_toggle").click()
    at.run()

    assert not at.exception
    assert "theme_target" not in at.session_state  # consumed by the script that applies it


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
