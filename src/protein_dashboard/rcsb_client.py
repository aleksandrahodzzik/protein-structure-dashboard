"""Thin client for the RCSB PDB REST APIs (Data API, file download, Search API v2)."""

from __future__ import annotations

import re
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import requests

DATA_API = "https://data.rcsb.org/rest/v1/core"
FILES_URL = "https://files.rcsb.org/download"
SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_PAGE = "https://www.rcsb.org/structure"

USER_AGENT = "protein-structure-dashboard/0.1 (+https://github.com/aleksandrahodzzik/protein-structure-dashboard)"
PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
MAX_RETRY_AFTER_SECONDS = 10.0
# Only the first few entities are looked up for the organism: large complexes
# (ribosomes) have dozens, and the header does not need all of them.
MAX_ORGANISM_ENTITIES = 3


class RcsbError(Exception):
    """Base class for all client errors."""


class InvalidPdbIdError(RcsbError, ValueError):
    """The string is not shaped like a PDB ID."""


class EntryNotFoundError(RcsbError):
    """RCSB has no entry with this ID."""


class RcsbUnavailableError(RcsbError):
    """RCSB could not be reached or kept failing after all retries."""


@dataclass(frozen=True)
class EntryMetadata:
    pdb_id: str
    title: str
    method: str
    resolution: float | None
    organisms: tuple[str, ...]
    deposition_date: str
    atom_count: int | None

    @property
    def url(self) -> str:
        return f"{ENTRY_PAGE}/{self.pdb_id}"


def is_valid_pdb_id(value: str) -> bool:
    return bool(PDB_ID_PATTERN.match(value.strip()))


def normalise_pdb_id(value: str) -> str:
    candidate = value.strip()
    if not PDB_ID_PATTERN.match(candidate):
        raise InvalidPdbIdError(
            f"'{value}' is not a valid PDB ID: expected 4 characters, a digit followed by "
            "three letters or digits (e.g. 1UBQ)."
        )
    return candidate.upper()


class RcsbClient:
    def __init__(
        self,
        session: requests.Session | None = None,
        *,
        timeout: tuple[float, float] = (5.0, 30.0),
        max_retries: int = 3,
        backoff_seconds: float = 0.5,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        self._session = session or requests.Session()
        self._session.headers["User-Agent"] = USER_AGENT
        self._timeout = timeout
        self._max_retries = max_retries
        self._backoff = backoff_seconds
        self._sleep = sleep

    def fetch_entry(self, pdb_id: str) -> EntryMetadata:
        """Entry-level metadata plus the organism(s) of the first polymer entities."""
        pdb_id = normalise_pdb_id(pdb_id)
        entry = self._get_json(f"{DATA_API}/entry/{pdb_id}", pdb_id)
        info = entry.get("rcsb_entry_info", {})
        resolutions = info.get("resolution_combined") or []
        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        return EntryMetadata(
            pdb_id=pdb_id,
            title=entry.get("struct", {}).get("title", "").strip(),
            method=(entry.get("exptl") or [{}])[0].get("method", "Unknown"),
            resolution=float(resolutions[0]) if resolutions else None,
            organisms=self._fetch_organisms(pdb_id, entity_ids[:MAX_ORGANISM_ENTITIES]),
            deposition_date=entry.get("rcsb_accession_info", {}).get("deposit_date", "")[:10],
            atom_count=info.get("deposited_atom_count"),
        )

    def fetch_structure(self, pdb_id: str) -> str:
        """The coordinate file as mmCIF text."""
        pdb_id = normalise_pdb_id(pdb_id)
        response = self._request("GET", f"{FILES_URL}/{pdb_id}.cif")
        self._raise_for_status(response, pdb_id)
        return response.text

    def search(self, query: str, limit: int = 10) -> list[str]:
        """Top entry IDs for a full-text keyword query, best match first."""
        query = query.strip()
        if not query:
            raise ValueError("Search query must not be empty.")
        payload = {
            "query": {"type": "terminal", "service": "full_text", "parameters": {"value": query}},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": limit}},
        }
        response = self._request("POST", SEARCH_API, json=payload)
        # The Search API answers "no hits" with 204 and an empty body.
        if response.status_code == 204:
            return []
        self._raise_for_status(response, query)
        return [hit["identifier"] for hit in response.json().get("result_set", [])]

    def _fetch_organisms(self, pdb_id: str, entity_ids: list[str]) -> tuple[str, ...]:
        names: list[str] = []
        for entity_id in entity_ids:
            try:
                entity = self._get_json(f"{DATA_API}/polymer_entity/{pdb_id}/{entity_id}", pdb_id)
            except RcsbError:
                # The header is still useful without an organism, so a failed
                # secondary lookup must not hide the whole entry.
                continue
            for source in entity.get("rcsb_entity_source_organism", []):
                name = source.get("ncbi_scientific_name")
                if name and name not in names:
                    names.append(name)
        return tuple(names)

    def _get_json(self, url: str, identifier: str) -> dict[str, Any]:
        response = self._request("GET", url)
        self._raise_for_status(response, identifier)
        data: dict[str, Any] = response.json()
        return data

    @staticmethod
    def _raise_for_status(response: requests.Response, identifier: str) -> None:
        if response.status_code == 404:
            raise EntryNotFoundError(f"RCSB has no entry for '{identifier}'.")
        if not response.ok:
            raise RcsbError(f"RCSB returned HTTP {response.status_code} for '{identifier}'.")

    def _request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        """Send a request, retrying timeouts, connection errors, 429 and 5xx with backoff."""
        reason = ""
        for attempt in range(self._max_retries + 1):
            delay = self._backoff * 2**attempt
            try:
                response = self._session.request(method, url, timeout=self._timeout, **kwargs)
            except (requests.Timeout, requests.ConnectionError) as exc:
                reason = type(exc).__name__
            else:
                if response.status_code not in RETRY_STATUSES:
                    return response
                reason = f"HTTP {response.status_code}"
                delay = max(delay, _retry_after_seconds(response))
            if attempt < self._max_retries:
                self._sleep(delay)
        raise RcsbUnavailableError(
            f"RCSB did not respond successfully after {self._max_retries + 1} attempts "
            f"(last error: {reason}). Please try again shortly."
        )


def _retry_after_seconds(response: requests.Response) -> float:
    try:
        return min(float(response.headers.get("Retry-After", 0)), MAX_RETRY_AFTER_SECONDS)
    except ValueError:  # Retry-After may also be an HTTP date; fall back to our backoff
        return 0.0
