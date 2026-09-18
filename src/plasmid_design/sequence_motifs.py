"""Heuristic cryptic-motif scanning, keyed by host domain (prokaryote vs
eukaryote, from each codon usage table's ``domain`` metadata field).

These are rule-based substring/simple-consensus screens, not an exhaustive
motif database or a trained splice-site predictor -- reported as such in
every result so the limitation is never silently hidden.
"""

from __future__ import annotations

from dataclasses import dataclass

# Internal Shine-Dalgarno-like ribosome-binding motifs (PRODUCTION mode:
# avoid these appearing internally, not just at the intended RBS).
SD_LIKE_MOTIFS: tuple[str, ...] = ("GGAGG", "AGGAGG")

# Eukaryotic cryptic regulatory/processing motifs.
POLYA_SIGNALS: tuple[str, ...] = ("AATAAA", "ATTAAA")
TATA_BOXES: tuple[str, ...] = ("TATAAA", "TATATAA")
AU_RICH_ELEMENT = "ATTTA"


@dataclass(frozen=True)
class MotifHit:
    motif: str
    kind: str  # e.g. "SD-like", "polyA signal", "TATA box", "AU-rich element", "splice donor/acceptor (heuristic)"
    start: int  # 1-indexed inclusive
    end: int  # 1-indexed inclusive


def _find_all(seq: str, motif: str, kind: str) -> list[MotifHit]:
    hits = []
    start = seq.find(motif)
    while start != -1:
        hits.append(MotifHit(motif=motif, kind=kind, start=start + 1, end=start + len(motif)))
        start = seq.find(motif, start + 1)
    return hits


def _find_splice_like(seq: str) -> list[MotifHit]:
    """Very simplified GT...AG consensus screen (splice donor 'GT' /
    acceptor 'AG' with minimal flanking): NOT an exhaustive splice-site
    predictor, just a heuristic warning for internal GT/AG dinucleotide
    pairs that resemble the canonical splice consensus."""

    hits = []
    n = len(seq)
    for i in range(n - 1):
        if seq[i : i + 2] == "GT":
            # look for a downstream AG within a plausible min-intron window
            window_end = min(n, i + 200)
            j = seq.find("AG", i + 20, window_end)
            if j != -1:
                hits.append(
                    MotifHit(
                        motif=seq[i : i + 2] + "..." + seq[j : j + 2],
                        kind="splice donor/acceptor (heuristic)",
                        start=i + 1,
                        end=j + 2,
                    )
                )
    return hits


def scan_motifs(dna: str, domain: str) -> list[MotifHit]:
    """Scan ``dna`` for domain-appropriate cryptic motifs. ``domain`` is
    'prokaryote' or 'eukaryote' (from CodonUsageTable.metadata['domain'])."""

    dna = dna.upper()
    hits: list[MotifHit] = []

    if domain == "prokaryote":
        for motif in SD_LIKE_MOTIFS:
            hits.extend(_find_all(dna, motif, "SD-like"))
    elif domain == "eukaryote":
        for motif in POLYA_SIGNALS:
            hits.extend(_find_all(dna, motif, "polyA signal"))
        for motif in TATA_BOXES:
            hits.extend(_find_all(dna, motif, "TATA box"))
        hits.extend(_find_all(dna, AU_RICH_ELEMENT, "AU-rich element"))
        hits.extend(_find_splice_like(dna))

    return hits
