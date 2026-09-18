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

from .dna_utils import SiteOccurrence, find_site_occurrences, reverse_complement

# (enzyme name, recognition sequence, severity)
RFC10_SITES: tuple[tuple[str, str, str], ...] = (
    ("EcoRI", "GAATTC", "fail"),
    ("XbaI", "TCTAGA", "fail"),
    ("SpeI", "ACTAGT", "fail"),
    ("PstI", "CTGCAG", "fail"),
    ("NotI", "GCGGCCGC", "warning"),
)

__all__ = [
    "RFC10_SITES",
    "reverse_complement",
    "find_site_occurrences",
    "SiteOccurrence",
    "SiteCheckResult",
    "RFC10Report",
    "validate_rfc10",
]


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
