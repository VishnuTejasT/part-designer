"""Forbidden restriction-site scanning beyond RFC10.

Golden Gate Type IIS enzymes are always screened (both strands), plus
RFC10's hard-fail sites, plus whatever the caller supplies for their own
cloning strategy. Built on the same IUPAC/both-strand scanner rfc10.py
uses (``dna_utils.find_site_occurrences``), so no scanning logic is
duplicated.
"""

from __future__ import annotations

import itertools

from .dna_utils import SiteOccurrence, find_site_occurrences
from .rfc10 import RFC10_SITES

# (enzyme name, recognition sequence) -- always screened regardless of the
# caller's own forbidden-enzyme list, since these are the standard Golden
# Gate / MoClo Type IIS enzymes most iGEM cloning strategies must avoid.
ALWAYS_ON_TYPE_IIS: tuple[tuple[str, str], ...] = (
    ("BsaI", "GGTCTC"),
    ("BsmBI", "CGTCTC"),
    ("BbsI", "GAAGAC"),
    ("SapI", "GCTCTTC"),
)


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
MAX_EXPANSIONS = 256


def expand_iupac(pattern: str) -> list[str]:
    """Concrete DNA sites matching a pattern that may use IUPAC ambiguity
    letters (R, Y, N, ...). Capped so a run of N's can't explode the scan."""

    try:
        choices = [IUPAC[base] for base in pattern.upper()]
    except KeyError as exc:
        raise ValueError(f"Invalid base {exc.args[0]!r} in cut site {pattern!r}") from None
    total = 1
    for c in choices:
        total *= len(c)
    if total > MAX_EXPANSIONS:
        raise ValueError(f"Cut site {pattern!r} is too ambiguous (matches {total} sequences; limit {MAX_EXPANSIONS})")
    return ["".join(combo) for combo in itertools.product(*choices)]


def scan_forbidden_sites(
    dna: str, extra_sites: tuple[tuple[str, str], ...] = ()
) -> dict[str, tuple[SiteOccurrence, ...]]:
    """Scan ``dna`` for RFC10 hard-fail sites, the always-on Type IIS
    enzymes, and any caller-supplied ``(name, pattern)`` pairs. Returns a
    dict of enzyme name -> occurrences (only enzymes with at least one hit
    are included)."""

    sites: list[tuple[str, str]] = [
        (name, pattern) for name, pattern, severity in RFC10_SITES if severity == "fail"
    ]
    sites.extend(ALWAYS_ON_TYPE_IIS)
    sites.extend(extra_sites)

    hits: dict[str, tuple[SiteOccurrence, ...]] = {}
    for name, pattern in sites:
        seen: set[tuple[int, int]] = set()
        for concrete in expand_iupac(pattern):
            for occ in find_site_occurrences(dna, concrete):
                key = (occ["start"], occ["end"])
                if key in seen:
                    continue
                seen.add(key)
                hits.setdefault(name, tuple())
                hits[name] = hits[name] + (SiteOccurrence(enzyme=name, pattern=pattern.upper(), **occ),)
    return hits
