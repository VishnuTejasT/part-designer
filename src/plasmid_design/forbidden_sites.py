"""Forbidden restriction-site scanning beyond RFC10.

Golden Gate Type IIS enzymes are always screened (both strands), plus
RFC10's hard-fail sites, plus whatever the caller supplies for their own
cloning strategy. Built on the same IUPAC/both-strand scanner rfc10.py
uses (``dna_utils.find_site_occurrences``), so no scanning logic is
duplicated.
"""

from __future__ import annotations

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
        raw = find_site_occurrences(dna, pattern)
        if raw:
            occurrences = tuple(
                SiteOccurrence(enzyme=name, pattern=pattern.upper(), **occ)
                for occ in raw
            )
            hits.setdefault(name, tuple())
            hits[name] = hits[name] + occurrences
    return hits
