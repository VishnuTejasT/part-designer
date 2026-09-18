"""BioBrick RFC10 assembly-standard compatibility validation.

Confirmed rule set (do not re-derive):

Illegal sites -- hard fail, checked on both strands:
    EcoRI  GAATTC
    XbaI   TCTAGA
    SpeI   ACTAGT
    PstI   CTGCAG

Flagged, not a hard fail:
    NotI   GCGGCCGC

Golden Gate / BsaI-BsmBI overhang validation is intentionally out of scope
for this check (a later addition), but the scanner below is written
generically (it does not assume a palindromic pattern) so a non-palindromic
site like BsaI's GGTCTC can be added to ``RFC10_SITES`` without changes to
the scanning logic.
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# (enzyme name, recognition sequence, severity)
RFC10_SITES: tuple[tuple[str, str, str], ...] = (
    ("EcoRI", "GAATTC", "fail"),
    ("XbaI", "TCTAGA", "fail"),
    ("SpeI", "ACTAGT", "fail"),
    ("PstI", "CTGCAG", "fail"),
    ("NotI", "GCGGCCGC", "warning"),
)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class SiteOccurrence:
    enzyme: str
    pattern: str
    strand: str  # "+", "-", or "+/-" (palindromic site, found on both scans)
    start: int  # 1-indexed, inclusive, forward-strand coordinates
    end: int  # 1-indexed, inclusive, forward-strand coordinates
    matched_sequence: str  # top-strand sequence spanning [start, end]


@dataclass(frozen=True)
class SiteCheckResult:
    enzyme: str
    pattern: str
    severity: str  # "fail" or "warning"
    occurrences: tuple[SiteOccurrence, ...]

    @property
    def passed(self) -> bool:
        return not self.occurrences


@dataclass(frozen=True)
class RFC10Report:
    dna: str
    checks: tuple[SiteCheckResult, ...]

    @property
    def hard_failures(self) -> tuple[SiteCheckResult, ...]:
        return tuple(
            c for c in self.checks if c.severity == "fail" and not c.passed
        )

    @property
    def warnings(self) -> tuple[SiteCheckResult, ...]:
        return tuple(
            c for c in self.checks if c.severity == "warning" and not c.passed
        )

    @property
    def passed(self) -> bool:
        """Overall pass/fail: warnings do not affect this, only hard fails."""
        return not self.hard_failures


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


def validate_rfc10(dna: str) -> RFC10Report:
    checks = []
    for enzyme, pattern, severity in RFC10_SITES:
        raw_occurrences = find_site_occurrences(dna, pattern)
        occurrences = tuple(
            SiteOccurrence(enzyme=enzyme, pattern=pattern, **occ)
            for occ in raw_occurrences
        )
        checks.append(
            SiteCheckResult(
                enzyme=enzyme,
                pattern=pattern,
                severity=severity,
                occurrences=occurrences,
            )
        )
    return RFC10Report(dna=dna.upper(), checks=tuple(checks))
