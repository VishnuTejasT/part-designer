"""tRNA Adaptation Index (tAI), dos Reis, Savva & Wernisch (2003, 2004).

Implements the exact algorithm and default wobble selective-constraint
(s) values from the canonical reference implementation (Mario dos Reis'
CRAN package ``tAI``, ``R/tAI.R``, function ``get.ws``/``get.tai``), verified
directly against the package source (CRAN ``tAI_0.2.2``) rather than
reconstructed from memory:

    s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)

References:
    dos Reis M., Wernisch L., Savva R. (2003) Unexpected correlations between
    gene expression and codon usage bias from microarray data for the whole
    Escherichia coli K-12 genome. Nucleic Acids Res. 31: 6976-85.
    dos Reis M., Savva R., Wernisch L. (2004) Solving the riddle of codon
    usage preferences: a test for translational selection. Nucleic Acids
    Res. 32: 5036-44.

Data requirement: tRNA gene copy numbers (GCN) per anticodon, sourced from
the Genomic tRNA Database (GtRNAdb, http://gtrnadb.ucsc.edu/) per host.
Only hosts with a real, sourced GCN table (see ``TRNA_TABLES``) can report
tAI; every other host must be reported as "tAI unavailable" -- this module
never fabricates a table.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache

from .codon_usage import DATA_DIR, CodonUsageTable

# Registry of hosts with a real, GtRNAdb-sourced tRNA gene-copy-number table.
# Same "drop in a JSON file, register it here" pattern as ORGANISM_TABLES in
# codon_usage.py -- extend this as more hosts get real tRNA data sourced.
TRNA_TABLES: dict[str, str] = {
    "e_coli_k12": "e_coli_k12_trna_gcn.json",
    "s_cerevisiae": "s_cerevisiae_trna_gcn.json",
    "b_subtilis_168": "b_subtilis_168_trna_gcn.json",
    "human": "human_trna_gcn.json",
}

# dos Reis et al. (2004) optimised selective-constraint (s) values, in the
# canonical order used by get.ws(): [WC x4, T:G, C:I, A:I, G:U, T:I]
DOS_REIS_S = (0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)

_BASES = "TCAG"
# Standard codonR/tAI codon ordering: nested T,C,A,G at each of the 3
# positions (index 0 = TTT, 1 = TTC, ..., 63 = GGG). This is the same order
# dos Reis' get.ws() assumes for its hardcoded box-of-4 wobble arithmetic.
CODON_ORDER: tuple[str, ...] = tuple(
    b1 + b2 + b3 for b1 in _BASES for b2 in _BASES for b3 in _BASES
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _anticodon_for(codon: str) -> str:
    """The Watson-Crick-complementary anticodon (5'->3') for a codon."""
    return codon.translate(_COMPLEMENT)[::-1]


class TrnaDataError(ValueError):
    pass


@dataclass(frozen=True)
class TrnaTable:
    organism: str
    metadata: dict
    anticodon_counts: dict[str, int]


@lru_cache(maxsize=None)
def load_trna_table(organism: str) -> TrnaTable:
    if organism not in TRNA_TABLES:
        available = ", ".join(sorted(TRNA_TABLES))
        raise TrnaDataError(
            f"No tRNA gene-copy-number table for {organism!r}. Available: "
            f"{available or '(none)'}"
        )
    path = DATA_DIR / TRNA_TABLES[organism]
    if not path.exists():
        raise TrnaDataError(f"tRNA GCN data file missing: {path}")
    with path.open() as fh:
        raw = json.load(fh)
    anticodons = {k: int(v) for k, v in raw["anticodons"].items()}
    metadata = {k: v for k, v in raw.items() if k != "anticodons"}
    return TrnaTable(organism=organism, metadata=metadata, anticodon_counts=anticodons)


def _gcn_vector(trna_table: TrnaTable) -> list[int]:
    """64-length tRNA GCN vector indexed like CODON_ORDER: entry i is the
    gene count of the anticodon exactly Watson-Crick-complementary to
    codon i (dos Reis' ``tRNA`` input convention)."""

    return [
        trna_table.anticodon_counts.get(_anticodon_for(codon), 0)
        for codon in CODON_ORDER
    ]


def relative_adaptiveness(
    trna_table: TrnaTable,
    codon_table: CodonUsageTable,
    s: tuple[float, ...] = DOS_REIS_S,
) -> dict[str, float]:
    """dos Reis et al.'s w (relative adaptiveness) per sense codon, excluding
    stops and Met (which get.ws() also drops). ``codon_table`` is used only
    to determine the host's super-kingdom (bacteria vs eukaryote) for the
    bacterial-only Ile-AUA lysidine special case, via its ``domain``
    metadata field."""

    p = [1.0 - si for si in s]
    tRNA = _gcn_vector(trna_table)
    domain = codon_table.metadata.get("domain")
    is_prokaryote = domain == "prokaryote"

    W = [0.0] * 64
    for i in range(0, 64, 4):
        # box of 4 codons sharing the first two bases, 3rd base T,C,A,G
        W[i] = p[0] * tRNA[i] + p[4] * tRNA[i + 1]  # ends T: WC + G:U wobble
        W[i + 1] = p[1] * tRNA[i + 1] + p[5] * tRNA[i]  # ends C: WC + I:C wobble
        W[i + 2] = p[2] * tRNA[i + 2] + p[6] * tRNA[i]  # ends A: WC + I:A wobble
        W[i + 3] = p[3] * tRNA[i + 3] + p[7] * tRNA[i + 2]  # ends G: WC + U:G wobble

    met_index = CODON_ORDER.index("ATG")
    W[met_index] = p[3] * tRNA[met_index]

    if is_prokaryote:
        ile_ata_index = CODON_ORDER.index("ATA")
        W[ile_ata_index] = p[8]  # bacterial lysidine-modified tRNA special case

    stop_codons = {"TAA", "TAG", "TGA"}
    kept = [
        (codon, w)
        for codon, w in zip(CODON_ORDER, W)
        if codon not in stop_codons and codon != "ATG"
    ]

    max_w = max(w for _, w in kept)
    if max_w <= 0:
        raise TrnaDataError(
            f"tRNA GCN table for {trna_table.organism!r} yields all-zero "
            "adaptiveness values -- cannot compute tAI"
        )
    normalized = {codon: w / max_w for codon, w in kept}

    zero_codons = [c for c, w in normalized.items() if w == 0]
    nonzero = [w for w in normalized.values() if w > 0]
    if zero_codons and nonzero:
        # Same imputation dos Reis' get.ws() uses: substitute the geometric
        # mean of the nonzero w's for codons with no decoding tRNA at all.
        geo_mean = math.exp(sum(math.log(w) for w in nonzero) / len(nonzero))
        for codon in zero_codons:
            normalized[codon] = geo_mean

    return normalized


def tai_score(dna: str, trna_table: TrnaTable, codon_table: CodonUsageTable) -> float:
    """Geometric mean of per-codon relative adaptiveness (w) across ``dna``,
    excluding stop codons and Met (per dos Reis' get.tai())."""

    w = relative_adaptiveness(trna_table, codon_table)
    dna = dna.upper()
    codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]

    log_sum = 0.0
    n = 0
    for codon in codons:
        if codon not in w:
            continue
        log_sum += math.log(w[codon])
        n += 1

    if n == 0:
        raise ValueError("No scorable codons (sequence is entirely stop/Met codons)")

    return math.exp(log_sum / n)
