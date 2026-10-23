import pytest
import requests
import responses

from protein_dashboard import rcsb_client
from protein_dashboard.rcsb_client import (
    DATA_API,
    FILES_URL,
    SEARCH_API,
    EntryNotFoundError,
    InvalidPdbIdError,
    RcsbClient,
    RcsbError,
    RcsbUnavailableError,
)

ENTRY_URL = f"{DATA_API}/entry/1UBQ"
ENTITY_URL = f"{DATA_API}/polymer_entity/1UBQ/1"
ENTRY_JSON = {
    "struct": {"title": "STRUCTURE OF UBIQUITIN REFINED AT 1.8 ANGSTROMS RESOLUTION "},
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "rcsb_entry_info": {"resolution_combined": [1.8], "deposited_atom_count": 660},
    "rcsb_accession_info": {"deposit_date": "1987-01-02T00:00:00.000+00:00"},
    "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"]},
}
ENTITY_JSON = {"rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens"}]}


@pytest.fixture
def sleeps() -> list[float]:
    return []


@pytest.fixture
def client(sleeps: list[float]) -> RcsbClient:
    return RcsbClient(max_retries=3, backoff_seconds=0.5, sleep=sleeps.append)


@pytest.mark.parametrize("value", ["1UBQ", "1ubq", " 4hhb ", "9ZZ9", "6VXX"])
def test_valid_ids_are_accepted_and_upper_cased(value: str) -> None:
    assert rcsb_client.normalise_pdb_id(value) == value.strip().upper()
    assert rcsb_client.is_valid_pdb_id(value)


@pytest.mark.parametrize("value", ["", "UBQ", "1UB", "1UBQ1", "AUBQ", "1UB!", "1 BQ", "ubiquitin"])
def test_invalid_ids_are_rejected(value: str) -> None:
    assert not rcsb_client.is_valid_pdb_id(value)
    with pytest.raises(InvalidPdbIdError, match="not a valid PDB ID"):
        rcsb_client.normalise_pdb_id(value)


@responses.activate
def test_fetch_entry_success_maps_fields(client: RcsbClient) -> None:
    responses.get(ENTRY_URL, json=ENTRY_JSON)
    responses.get(ENTITY_URL, json=ENTITY_JSON)

    entry = client.fetch_entry("1ubq")

    assert entry.pdb_id == "1UBQ"
    assert entry.title == "STRUCTURE OF UBIQUITIN REFINED AT 1.8 ANGSTROMS RESOLUTION"
    assert entry.method == "X-RAY DIFFRACTION"
    assert entry.resolution == 1.8
    assert entry.organisms == ("Homo sapiens",)
    assert entry.deposition_date == "1987-01-02"
    assert entry.atom_count == 660
    assert entry.url == "https://www.rcsb.org/structure/1UBQ"


@responses.activate
def test_fetch_entry_can_skip_organism_lookups(client: RcsbClient) -> None:
    responses.get(ENTRY_URL, json=ENTRY_JSON)

    entry = client.fetch_entry("1UBQ", with_organisms=False)

    assert entry.organisms == ()
    assert len(responses.calls) == 1


@responses.activate
def test_requests_send_a_user_agent_and_timeout(client: RcsbClient) -> None:
    responses.get(ENTRY_URL, json=ENTRY_JSON)
    responses.get(ENTITY_URL, json=ENTITY_JSON)

    client.fetch_entry("1UBQ")

    request = responses.calls[0].request
    assert request.headers["User-Agent"].startswith("protein-structure-dashboard/")
    assert request.req_kwargs["timeout"] == (5.0, 30.0)


@responses.activate
def test_fetch_entry_without_resolution_or_organism(client: RcsbClient) -> None:
    nmr = {**ENTRY_JSON, "rcsb_entry_info": {}, "exptl": [{"method": "SOLUTION NMR"}]}
    responses.get(ENTRY_URL, json=nmr)
    responses.get(ENTITY_URL, status=404)

    entry = client.fetch_entry("1UBQ")

    assert entry.resolution is None
    assert entry.atom_count is None
    assert entry.organisms == ()


@responses.activate
def test_fetch_entry_unknown_id_raises_not_found(client: RcsbClient) -> None:
    responses.get(f"{DATA_API}/entry/0ZZZ", status=404, json={"status": 404})

    with pytest.raises(EntryNotFoundError, match="0ZZZ"):
        client.fetch_entry("0ZZZ")
    assert len(responses.calls) == 1  # a 404 is final and must not be retried


def test_fetch_entry_invalid_id_makes_no_request(client: RcsbClient) -> None:
    with pytest.raises(InvalidPdbIdError):
        client.fetch_entry("nope")


@responses.activate
def test_fetch_structure_returns_mmcif_text(client: RcsbClient) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", body="data_1UBQ\n#\n")

    assert client.fetch_structure("1UBQ").startswith("data_1UBQ")


@responses.activate
def test_fetch_structure_404(client: RcsbClient) -> None:
    responses.get(f"{FILES_URL}/0ZZZ.cif", status=404, body="<html>Not Found</html>")

    with pytest.raises(EntryNotFoundError):
        client.fetch_structure("0ZZZ")


@responses.activate
def test_429_then_success_is_retried_with_backoff(client: RcsbClient, sleeps: list[float]) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", status=429)
    responses.get(f"{FILES_URL}/1UBQ.cif", status=503)
    responses.get(f"{FILES_URL}/1UBQ.cif", body="data_1UBQ")

    assert client.fetch_structure("1UBQ") == "data_1UBQ"
    assert len(responses.calls) == 3
    assert sleeps == [0.5, 1.0]  # exponential backoff


@responses.activate
def test_retry_after_header_is_honoured_but_capped(client: RcsbClient, sleeps: list[float]) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", status=429, headers={"Retry-After": "3"})
    responses.get(f"{FILES_URL}/1UBQ.cif", status=429, headers={"Retry-After": "999"})
    responses.get(f"{FILES_URL}/1UBQ.cif", body="data_1UBQ")

    client.fetch_structure("1UBQ")

    assert sleeps == [3.0, rcsb_client.MAX_RETRY_AFTER_SECONDS]


@responses.activate
def test_persistent_server_error_gives_up_with_clear_message(
    client: RcsbClient, sleeps: list[float]
) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", status=500)

    with pytest.raises(RcsbUnavailableError, match=r"4 attempts.*HTTP 500"):
        client.fetch_structure("1UBQ")
    assert len(responses.calls) == 4
    assert len(sleeps) == 3  # no sleep after the final attempt


@responses.activate
def test_timeout_is_retried_then_reported(client: RcsbClient) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", body=requests.exceptions.ConnectTimeout())

    with pytest.raises(RcsbUnavailableError, match="ConnectTimeout"):
        client.fetch_structure("1UBQ")
    assert len(responses.calls) == 4


@responses.activate
def test_timeout_then_success(client: RcsbClient) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", body=requests.exceptions.ReadTimeout())
    responses.get(f"{FILES_URL}/1UBQ.cif", body="data_1UBQ")

    assert client.fetch_structure("1UBQ") == "data_1UBQ"


@responses.activate
def test_other_client_errors_are_not_retried(client: RcsbClient) -> None:
    responses.get(f"{FILES_URL}/1UBQ.cif", status=403)

    with pytest.raises(RcsbError, match="HTTP 403"):
        client.fetch_structure("1UBQ")
    assert len(responses.calls) == 1


@responses.activate
def test_search_returns_identifiers_in_rank_order(client: RcsbClient) -> None:
    hits = [{"identifier": "1UBQ", "score": 1.0}, {"identifier": "1UBI", "score": 0.9}]
    responses.post(SEARCH_API, json={"total_count": 2, "result_set": hits})

    assert client.search(" ubiquitin ", limit=10) == ["1UBQ", "1UBI"]

    body = responses.calls[0].request.body
    assert body is not None
    assert b'"value": "ubiquitin"' in body
    assert b'"rows": 10' in body


@responses.activate
def test_search_with_no_hits_returns_empty_list(client: RcsbClient) -> None:
    responses.post(SEARCH_API, status=204)

    assert client.search("zzzqqqxxx") == []


def test_search_rejects_blank_query(client: RcsbClient) -> None:
    with pytest.raises(ValueError, match="empty"):
        client.search("   ")
