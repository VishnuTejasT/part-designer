"""Codon Adaptation Index (CAI).

Per-codon weight w = usage(codon) / usage(most-used synonymous codon) for
that amino acid. Since ``CodonOption.fraction`` is already usage relative to
the sum of synonymous codons for the same amino acid, w = fraction /
max(fraction) within the amino acid group is mathematically identical to
using raw per-thousand usage for the ratio (both are proportional within a
group), so no new data is needed beyond the existing codon usage tables.

CAI is the geometric mean of w across the sequence, excluding Met, Trp
(both have a single codon, so w is trivially 1.0 and they carry no
information) and the stop codon.
"""

from __future__ import annotations

import math

from .codon_usage import CodonUsageTable

EXCLUDED_AMINO_ACIDS = {"M", "W", "*"}


def w_scores(table: CodonUsageTable) -> dict[str, float]:
    """Map every codon in ``table`` to its w = fraction / max_fraction
    within its synonymous group."""

    scores: dict[str, float] = {}
    for amino_acid, options in table.by_amino_acid.items():
        max_fraction = max(opt.fraction for opt in options)
        for opt in options:
            scores[opt.codon] = (
                opt.fraction / max_fraction if max_fraction > 0 else 0.0
            )
    return scores


def cai_score(dna: str, table: CodonUsageTable) -> float:
    """Geometric mean of per-codon w across ``dna``, excluding Met/Trp/stop.

    Raises ValueError if no scorable codons remain (e.g. a 1-2 residue
    protein consisting solely of Met/Trp/stop).
    """

    scores = w_scores(table)
    dna = dna.upper()
    codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]

    log_sum = 0.0
    n = 0
    for codon in codons:
        info = table.by_codon.get(codon)
        if info is None:
            continue
        if info["amino_acid"] in EXCLUDED_AMINO_ACIDS:
            continue
        w = scores.get(codon, 0.0)
        if w <= 0:
            # A codon with essentially zero usage in this host is a strong
            # penalty, not an undefined log -- floor it rather than blow up.
            w = 1e-6
        log_sum += math.log(w)
        n += 1

    if n == 0:
        raise ValueError(
            "No scorable codons (sequence consists solely of Met/Trp/stop "
            "codons, which are excluded from the CAI definition)"
        )

    return math.exp(log_sum / n)


def expression_weighted_cai_score(dna: str, table: CodonUsageTable) -> tuple[None, str]:
    """No expression-weighted (highly-expressed-gene-subset) codon usage
    table has been sourced for any host yet -- always returns (None, reason)
    rather than silently omitting the field or faking a number."""

    return None, (
        f"no expression-weighted codon usage table sourced for "
        f"{table.organism!r}; only the whole-genome table is available"
    )
