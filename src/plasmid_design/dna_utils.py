"""Sequence-level DNA primitives shared across restriction-site scanning,
forbidden-site scanning, and hard-constraint checks (GC windows, homopolymer
runs, repeated k-mers). Nothing here is organism- or enzyme-specific.
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base in "GC")
    return gc / len(seq)


@dataclass(frozen=True)
class SiteOccurrence:
    enzyme: str
    pattern: str
    strand: str  # "+", "-", or "+/-" (palindromic site, found on both scans)
    start: int  # 1-indexed, inclusive, forward-strand coordinates
    end: int  # 1-indexed, inclusive, forward-strand coordinates
    matched_sequence: str  # top-strand sequence spanning [start, end]


def _forward_matches(seq: str, pattern: str) -> set[tuple[int, int]]:
    """1-indexed inclusive (start, end) spans of every (overlapping)
    occurrence of pattern in seq."""

    k = len(pattern)
    spans = set()
    start = seq.find(pattern)
    while start != -1:
        spans.add((start + 1, start + k))
        start = seq.find(pattern, start + 1)
    return spans


def find_site_occurrences(seq: str, pattern: str) -> list[dict]:
    """Scan both strands of ``seq`` for ``pattern``, returning forward-strand
    coordinates for every hit (deduplicated when a palindromic pattern is
    found at the same location from both scans)."""

    seq = seq.upper()
    pattern = pattern.upper()
    k = len(pattern)
    length = len(seq)

    forward_spans = _forward_matches(seq, pattern)

    rc_seq = reverse_complement(seq)
    reverse_spans = set()
    for i, _ in enumerate(rc_seq):
        if rc_seq[i : i + k] == pattern:
            start_1 = length - i - k + 1
            end_1 = length - i
            reverse_spans.add((start_1, end_1))

    all_spans = forward_spans | reverse_spans
    occurrences = []
    for start, end in sorted(all_spans):
        on_forward = (start, end) in forward_spans
        on_reverse = (start, end) in reverse_spans
        if on_forward and on_reverse:
            strand = "+/-"
        elif on_forward:
            strand = "+"
        else:
            strand = "-"
        occurrences.append(
            {
                "strand": strand,
                "start": start,
                "end": end,
                "matched_sequence": seq[start - 1 : end],
            }
        )
    return occurrences


@dataclass(frozen=True)
class HomopolymerRun:
    base: str
    start: int  # 1-indexed inclusive
    end: int  # 1-indexed inclusive
    length: int


def homopolymer_runs(seq: str, max_len: int = 6) -> list[HomopolymerRun]:
    """Runs of a single repeated base longer than ``max_len``."""

    seq = seq.upper()
    runs = []
    i = 0
    n = len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        run_len = j - i
        if run_len > max_len:
            runs.append(
                HomopolymerRun(base=seq[i], start=i + 1, end=j, length=run_len)
            )
        i = j
    return runs


@dataclass(frozen=True)
class RepeatedKmer:
    kmer: str
    positions: tuple[int, ...]  # 1-indexed starts of every occurrence


def repeated_kmers(seq: str, k: int = 15) -> list[RepeatedKmer]:
    """k-mers (default 15) that occur more than once in ``seq``."""

    seq = seq.upper()
    positions: dict[str, list[int]] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        positions.setdefault(kmer, []).append(i + 1)
    return [
        RepeatedKmer(kmer=kmer, positions=tuple(starts))
        for kmer, starts in positions.items()
        if len(starts) > 1
    ]


@dataclass(frozen=True)
class GcWindowViolation:
    start: int  # 1-indexed inclusive
    end: int  # 1-indexed inclusive
    gc_fraction: float
    bound: str  # "low" or "high"


def gc_window_scan(
    seq: str, window: int = 50, low: float = 0.30, high: float = 0.70
) -> list[GcWindowViolation]:
    """Sliding-window GC% scan; returns every window outside [low, high]."""

    seq = seq.upper()
    n = len(seq)
    if n < window:
        frac = gc_content(seq)
        if frac < low:
            return [GcWindowViolation(1, n, frac, "low")]
        if frac > high:
            return [GcWindowViolation(1, n, frac, "high")]
        return []

    violations = []
    for i in range(0, n - window + 1):
        chunk = seq[i : i + window]
        frac = gc_content(chunk)
        if frac < low:
            violations.append(GcWindowViolation(i + 1, i + window, frac, "low"))
        elif frac > high:
            violations.append(GcWindowViolation(i + 1, i + window, frac, "high"))
    return violations
