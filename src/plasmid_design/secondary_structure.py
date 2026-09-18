"""RNA secondary-structure hard-constraint checks (hairpins/stems/
terminator-like motifs/inverted repeats) via ViennaRNA, at a caller-supplied
folding temperature. Lazy ``import RNA`` (see folding.py's rationale) --
every function here returns ``None``/``[]`` rather than raising if
ViennaRNA is unavailable, and the caller records why.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dna_utils import gc_content, reverse_complement

DEFAULT_TEMPERATURE_C = 37.0
WINDOW = 60
STEP = 6
MAX_STEM_BP = 8  # hard constraint 6b: no stem >= this many contiguous bp
MAX_WINDOW_MFE = -15.0  # hard constraint 6c: no window more stable than this
INIT_DG_FLOOR = -5.0  # hard constraint 6a


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _fold(seq: str, temperature_c: float) -> tuple[str, float] | None:
    try:
        import RNA  # noqa: PLC0415
    except ImportError:
        return None
    try:
        RNA.cvar.temperature = float(temperature_c)
        structure, energy = RNA.fold(_dna_to_rna(seq))
    except Exception:
        return None
    return structure, float(energy)


def fold(seq: str, temperature_c: float = DEFAULT_TEMPERATURE_C) -> tuple[str, float] | None:
    """Public wrapper around the lazy ViennaRNA fold, for callers (like the
    optimizer's incremental refold) that need to fold an arbitrary slice."""

    return _fold(seq, temperature_c)


def _ptable(structure: str) -> list[int] | None:
    try:
        import RNA  # noqa: PLC0415
    except ImportError:
        return None
    return list(RNA.ptable(structure))


@dataclass(frozen=True)
class WindowFold:
    start: int  # 1-indexed inclusive, DNA coordinates
    end: int
    structure: str
    mfe: float


def sliding_window_fold(
    dna: str,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    window: int = WINDOW,
    step: int = STEP,
) -> list[WindowFold] | None:
    """Fold every ``window``-nt slice of ``dna`` at ``step``-nt increments.
    Returns None if ViennaRNA is unavailable."""

    n = len(dna)
    if n == 0:
        return []
    results = []
    starts = list(range(0, max(1, n - window + 1), step))
    if not starts or starts[-1] + window < n:
        starts.append(max(0, n - window))
    seen = set()
    for start in starts:
        if start in seen:
            continue
        seen.add(start)
        chunk = dna[start : start + window]
        if len(chunk) < 2:
            continue
        folded = _fold(chunk, temperature_c)
        if folded is None:
            return None
        structure, mfe = folded
        results.append(WindowFold(start=start + 1, end=start + len(chunk), structure=structure, mfe=mfe))
    return results


def longest_stem_length(structure: str, ptable: list[int] | None = None) -> int:
    """Longest run of contiguous base pairs (a single unbroken helix) in a
    dot-bracket structure -- NOT the same as longest run of '(' characters,
    which can span a bulge/interior loop."""

    pt = ptable if ptable is not None else _ptable(structure)
    if pt is None:
        return 0
    n = pt[0]
    best = 0
    i = 1
    while i <= n:
        if pt[i] > i:
            run = 1
            j = i
            while j + 1 <= n and pt[j + 1] == pt[j] - 1 and pt[j + 1] > j + 1:
                run += 1
                j += 1
            best = max(best, run)
            i = j + 1
        else:
            i += 1
    return best


def worst_window_metrics(windows: list[WindowFold]) -> tuple[float | None, int]:
    """(most negative MFE, longest stem) across all windows."""

    if not windows:
        return None, 0
    worst_mfe = min(w.mfe for w in windows)
    longest = max(longest_stem_length(w.structure) for w in windows)
    return worst_mfe, longest


@dataclass(frozen=True)
class TerminatorMotif:
    window_start: int
    stem_start: int  # 1-indexed, DNA coordinates
    stem_end: int
    stem_length: int
    stem_gc: float
    t_run_start: int
    t_run_length: int


def terminator_like_motifs(dna: str, windows: list[WindowFold]) -> list[TerminatorMotif]:
    """A GC-rich stem of 6+ bp followed by a run of 4+ T within 15 nt of the
    stem's end -- the classic intrinsic (rho-independent) terminator shape,
    screened as a hard constraint in prokaryotic hosts."""

    hits: list[TerminatorMotif] = []
    dna = dna.upper()
    for w in windows:
        pt = _ptable(w.structure)
        if pt is None:
            continue
        n = pt[0]
        i = 1
        while i <= n:
            if pt[i] > i:
                run_start = i
                j = i
                while j + 1 <= n and pt[j + 1] == pt[j] - 1 and pt[j + 1] > j + 1:
                    j += 1
                stem_len = j - run_start + 1
                if stem_len >= 6:
                    # The opening (5') arm runs run_start..j; its partners
                    # pt[run_start]..pt[j] are the closing (3') arm, with
                    # pt[run_start] the outermost/3'-most base of the whole
                    # hairpin -- that is where an unpaired poly-U tail would
                    # sit in vivo, not after the opening arm itself.
                    abs_stem_start = w.start + run_start - 1
                    abs_hairpin_end = w.start + pt[run_start] - 1
                    opening_seq = dna[abs_stem_start - 1 : w.start + j - 1]
                    closing_seq = dna[w.start + pt[j] - 1 - 1 : abs_hairpin_end]
                    stem_seq = opening_seq + closing_seq
                    if gc_content(stem_seq) >= 0.6:
                        downstream = dna[abs_hairpin_end : abs_hairpin_end + 15]
                        t_run = 0
                        best_run = 0
                        best_start = -1
                        for k, base in enumerate(downstream):
                            if base == "T":
                                t_run += 1
                                if t_run > best_run:
                                    best_run = t_run
                                    best_start = abs_hairpin_end + 1 + k - t_run + 1
                            else:
                                t_run = 0
                        if best_run >= 4:
                            hits.append(
                                TerminatorMotif(
                                    window_start=w.start,
                                    stem_start=abs_stem_start,
                                    stem_end=abs_hairpin_end,
                                    stem_length=stem_len,
                                    stem_gc=gc_content(stem_seq),
                                    t_run_start=best_start,
                                    t_run_length=best_run,
                                )
                            )
                i = j + 1
            else:
                i += 1
    return hits


@dataclass(frozen=True)
class InvertedRepeat:
    left_start: int  # 1-indexed
    left_end: int
    right_start: int
    right_end: int
    length: int


def inverted_repeats(dna: str, min_len: int = 8, max_span: int = 30) -> list[InvertedRepeat]:
    """Perfect inverted repeats of >= min_len bp whose two arms both fall
    within a max_span-nt window (i.e. could fold back into a hairpin)."""

    dna = dna.upper()
    n = len(dna)
    hits = []
    for i in range(n - min_len + 1):
        left = dna[i : i + min_len]
        rc = reverse_complement(left)
        # right arm must start after the left arm ends, and its start must
        # fall within max_span nt of the left arm's start (str.find's `end`
        # bounds where a match must fully fit, not just where it can start,
        # so the span check is done manually after each find).
        search_start = i + min_len
        pos = dna.find(rc, search_start)
        while pos != -1 and pos < i + max_span:
            hits.append(
                InvertedRepeat(
                    left_start=i + 1, left_end=i + min_len,
                    right_start=pos + 1, right_end=pos + min_len,
                    length=min_len,
                )
            )
            pos = dna.find(rc, pos + 1)
    return hits
