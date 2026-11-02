"""Sequence analytics. Pure functions: no I/O, no Streamlit, trivially unit-testable."""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable
from dataclasses import dataclass

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

AMINO_ACID_NAMES: dict[str, str] = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartate",
    "E": "Glutamate",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}

PROPERTY_CLASSES: dict[str, str] = {
    "Hydrophobic": "AVLIMFW",
    "Polar": "STNQY",
    "Charged": "DEKRH",
    "Special": "CGP",
}

# Window 19 with a 1.6 cut-off is the classic Kyte-Doolittle transmembrane-helix heuristic.
MIN_WINDOW = 3
MAX_WINDOW = 39
DEFAULT_WINDOW = 9
TM_WINDOW = 19
TM_THRESHOLD = 1.6


@dataclass(frozen=True)
class CleanedSequence:
    """A sequence reduced to the 20 standard amino acids, plus what was dropped."""

    sequence: str
    unknown_count: int


@dataclass(frozen=True)
class HydrophobicityProfile:
    """Sliding-window Kyte-Doolittle scores.

    Positions are the 1-based window centres, as indices into the analysed sequence.
    Within a gap-free stretch they are consecutive, so a jump of more than one between
    neighbouring positions marks a gap that no window was allowed to cross.
    """

    window: int
    positions: list[int]
    values: list[float]


@dataclass(frozen=True)
class ProteinMetrics:
    length: int
    unknown_count: int
    molecular_weight: float
    isoelectric_point: float
    gravy: float
    instability_index: float
    aromaticity: float
    composition: dict[str, float]


def _normalise(sequence: str) -> str:
    return "".join(sequence.split()).upper()


def clean_sequence(sequence: str) -> CleanedSequence:
    """Upper-case, drop whitespace, and strip anything that is not a standard residue.

    Biopython raises on X, B, Z, U, ... so non-standard letters are excluded
    from every calculation and reported as a count instead of silently guessed.
    """
    upper = _normalise(sequence)
    kept = "".join(c for c in upper if c in STANDARD_AMINO_ACIDS)
    return CleanedSequence(sequence=kept, unknown_count=len(upper) - len(kept))


def _analyser(
    sequence: str, *, monoisotopic: bool = False
) -> tuple[ProteinAnalysis, CleanedSequence]:
    cleaned = clean_sequence(sequence)
    if not cleaned.sequence:
        raise ValueError("Sequence has no standard amino acid residues to analyse.")
    return ProteinAnalysis(cleaned.sequence, monoisotopic=monoisotopic), cleaned


def composition_percent(sequence: str) -> dict[str, float]:
    """Percentage of each of the 20 standard residues (sums to 100)."""
    analyser, _ = _analyser(sequence)
    percent = analyser.amino_acids_percent
    return {aa: percent[aa] for aa in STANDARD_AMINO_ACIDS}


def molecular_weight(sequence: str, *, monoisotopic: bool = False) -> float:
    """Molecular weight in Da. Average masses by default, as reported by most tools."""
    analyser, _ = _analyser(sequence, monoisotopic=monoisotopic)
    return float(analyser.molecular_weight())


def contiguous_runs(sequence: str, breaks: Iterable[int] = ()) -> list[tuple[int, int]]:
    """Half-open ``(start, end)`` index ranges of consecutive standard residues.

    A run ends at any non-standard letter (X, B, ...) and after every index ``i`` in
    ``breaks``, which marks a gap between residue ``i`` and ``i + 1`` (for example
    residues missing from the electron density). Indices refer to the sequence with
    whitespace removed.
    """
    letters = _normalise(sequence)
    cuts = set(breaks)
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, letter in enumerate(letters):
        if letter not in STANDARD_AMINO_ACIDS:
            if start is not None:
                runs.append((start, index))
            start = None
            continue
        if start is None:
            start = index
        if index in cuts:
            runs.append((start, index + 1))
            start = None
    if start is not None:
        runs.append((start, len(letters)))
    return runs


def longest_run(sequence: str, breaks: Iterable[int] = ()) -> int:
    """Length of the longest gap-free stretch of standard residues, 0 if there is none."""
    return max((end - start for start, end in contiguous_runs(sequence, breaks)), default=0)


def kd_profile(
    sequence: str, window: int = DEFAULT_WINDOW, *, breaks: Iterable[int] = ()
) -> HydrophobicityProfile:
    """Kyte-Doolittle hydropathy averaged over a sliding window.

    The window is centred on each residue, so it must be odd. A window is only scored
    when all of its residues are standard and contiguous: it never crosses a gap listed
    in ``breaks`` or a non-standard residue, because averaging residues that are not
    neighbours in the chain would invent a hydrophobic or polar stretch. Each gap-free
    stretch of length ``n`` contributes ``n - window + 1`` points, and stretches shorter
    than the window contribute none.
    """
    letters = _normalise(sequence)
    if not clean_sequence(letters).sequence:
        raise ValueError("Sequence has no standard amino acid residues to analyse.")
    if window < 3 or window % 2 == 0:
        raise ValueError(f"Window must be an odd integer >= 3, got {window}.")
    runs = contiguous_runs(letters, breaks)
    longest = max(end - start for start, end in runs)
    if window > longest:
        raise ValueError(
            f"Window {window} is longer than the longest gap-free stretch ({longest} residues)."
        )
    half = window // 2
    positions: list[int] = []
    values: list[float] = []
    for start, end in runs:
        if end - start < window:
            continue
        scores = ProteinAnalysis(letters[start:end]).protein_scale(kd, window)
        positions.extend(start + half + 1 + i for i in range(len(scores)))
        values.extend(float(v) for v in scores)
    return HydrophobicityProfile(window=window, positions=positions, values=values)


def max_window(length: int) -> int:
    """Largest valid (odd) window for a sequence of ``length`` residues, 0 if none fits."""
    limit = min(length, MAX_WINDOW)
    if limit % 2 == 0:
        limit -= 1
    return limit if limit >= MIN_WINDOW else 0


def compute_metrics(sequence: str) -> ProteinMetrics:
    analyser, cleaned = _analyser(sequence)
    percent = analyser.amino_acids_percent
    return ProteinMetrics(
        length=len(cleaned.sequence),
        unknown_count=cleaned.unknown_count,
        molecular_weight=float(analyser.molecular_weight()),
        isoelectric_point=float(analyser.isoelectric_point()),
        gravy=float(analyser.gravy()),
        instability_index=float(analyser.instability_index()),
        aromaticity=float(analyser.aromaticity()),
        composition={aa: percent[aa] for aa in STANDARD_AMINO_ACIDS},
    )


def residue_kd_value(one_letter: str) -> float | None:
    """Per-residue Kyte-Doolittle value, or None for non-standard residues."""
    return kd.get(one_letter.upper())


def metrics_to_csv(label: str, metrics: ProteinMetrics) -> str:
    """Flat two-section CSV (summary metrics, then composition) for download."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["section", "name", "value"])
    summary = {
        "length_residues": metrics.length,
        "unknown_residues_excluded": metrics.unknown_count,
        "molecular_weight_da": round(metrics.molecular_weight, 2),
        "isoelectric_point": round(metrics.isoelectric_point, 2),
        "gravy": round(metrics.gravy, 3),
        "instability_index": round(metrics.instability_index, 2),
        "aromaticity": round(metrics.aromaticity, 3),
    }
    for name, value in summary.items():
        writer.writerow(["summary", f"{label}:{name}", value])
    for aa, pct in metrics.composition.items():
        writer.writerow(["composition_percent", f"{label}:{aa}", round(pct, 2)])
    return out.getvalue()
