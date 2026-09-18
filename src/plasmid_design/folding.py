"""mRNA secondary-structure folding via ViennaRNA (real MFE, not a
heuristic). ``RNA`` is imported lazily inside each function -- if the
ViennaRNA Python bindings aren't installed or fail to load in some
environment, every function here returns ``None`` rather than raising, and
the caller (optimize.py) records *why* in the result's notes. This keeps
the rest of the package usable even where ViennaRNA can't run.
"""

from __future__ import annotations


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def mfe(seq: str) -> tuple[str, float] | None:
    """Minimum free energy secondary structure and its dG (kcal/mol) for a
    DNA or RNA sequence. Returns None if ViennaRNA is unavailable."""

    try:
        import RNA  # noqa: PLC0415 -- intentionally lazy, see module docstring
    except ImportError:
        return None

    try:
        structure, energy = RNA.fold(_dna_to_rna(seq))
    except Exception:
        return None
    return structure, float(energy)


def initiation_window(
    cds: str, five_prime_utr: str, upstream: int = 20, downstream: int = 40
) -> str:
    """The -upstream..+downstream window around the start codon, per the
    spec: 5' UTR (if shorter than ``upstream``, padded with what's
    available) plus the first ``downstream`` bases of the CDS."""

    utr_slice = five_prime_utr[-upstream:] if five_prime_utr else ""
    cds_slice = cds[:downstream]
    return utr_slice + cds_slice


def initiation_dG(cds: str, five_prime_utr: str) -> float | None:
    """dG of the -20/+40 initiation window. Higher (less negative) dG means
    weaker secondary structure, which is what PRODUCTION mode wants to
    maximize. Returns None if ViennaRNA is unavailable."""

    window = initiation_window(cds, five_prime_utr)
    if not window:
        return None
    result = mfe(window)
    if result is None:
        return None
    return result[1]
