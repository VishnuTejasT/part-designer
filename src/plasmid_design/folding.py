"""mRNA secondary-structure folding via ViennaRNA (real MFE, not a
heuristic). ``RNA`` is imported lazily inside each function -- if the
ViennaRNA Python bindings aren't installed or fail to load in some
environment, every function here returns ``None`` rather than raising, and
the caller (optimize.py) records *why* in the result's notes. This keeps
the rest of the package usable even where ViennaRNA can't run.
"""

from __future__ import annotations

DEFAULT_TEMPERATURE_C = 37.0
UPSTREAM = 20
DOWNSTREAM = 40
UNPAIRED_UPSTREAM = 15  # hard constraint 6a: -15..+20 must be single-stranded
UNPAIRED_DOWNSTREAM = 20


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def mfe(seq: str, temperature_c: float = DEFAULT_TEMPERATURE_C) -> tuple[str, float] | None:
    """Minimum free energy secondary structure and its dG (kcal/mol) for a
    DNA or RNA sequence, at the given folding temperature (Celsius).
    Returns None if ViennaRNA is unavailable."""

    try:
        import RNA  # noqa: PLC0415 -- intentionally lazy, see module docstring
    except ImportError:
        return None

    try:
        RNA.cvar.temperature = float(temperature_c)
        structure, energy = RNA.fold(_dna_to_rna(seq))
    except Exception:
        return None
    return structure, float(energy)


def initiation_window(
    cds: str, five_prime_utr: str, upstream: int = UPSTREAM, downstream: int = DOWNSTREAM
) -> str:
    """The -upstream..+downstream window around the start codon, per the
    spec: 5' UTR (if shorter than ``upstream``, padded with what's
    available) plus the first ``downstream`` bases of the CDS."""

    utr_slice = five_prime_utr[-upstream:] if five_prime_utr else ""
    cds_slice = cds[:downstream]
    return utr_slice + cds_slice


def initiation_dG(
    cds: str, five_prime_utr: str, temperature_c: float = DEFAULT_TEMPERATURE_C
) -> float | None:
    """dG of the -20/+40 initiation window. Higher (less negative) dG means
    weaker secondary structure, which is what PRODUCTION mode wants to
    maximize toward ~0. Returns None if ViennaRNA is unavailable."""

    window = initiation_window(cds, five_prime_utr)
    if not window:
        return None
    result = mfe(window, temperature_c)
    if result is None:
        return None
    return result[1]


def initiation_unpaired_ok(
    cds: str,
    five_prime_utr: str,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    upstream: int = UPSTREAM,
    downstream: int = DOWNSTREAM,
) -> bool | None:
    """Hard constraint 6a: within the -20/+40 window, the RBS/Kozak and
    start codon region (-15..+20 relative to the start codon) must be fully
    single-stranded in the MFE structure. Returns None if ViennaRNA is
    unavailable or the UTR is too short to evaluate."""

    window = initiation_window(cds, five_prime_utr, upstream, downstream)
    if not window:
        return None
    result = mfe(window, temperature_c)
    if result is None:
        return None
    structure, _energy = result

    try:
        import RNA  # noqa: PLC0415

        pt = list(RNA.ptable(structure))
    except ImportError:
        return None

    # The window is utr_slice (up to `upstream` nt) + cds[:downstream]; the
    # start codon sits right after the UTR slice. Map the -15..+20 region
    # into window-local (1-indexed) coordinates.
    utr_len = min(len(five_prime_utr), upstream)
    start_codon_pos = utr_len  # 0-indexed offset of the start codon in window
    region_start = max(1, start_codon_pos - UNPAIRED_UPSTREAM + 1)
    region_end = min(len(window), start_codon_pos + UNPAIRED_DOWNSTREAM)

    return all(pt[i] == 0 for i in range(region_start, region_end + 1))
