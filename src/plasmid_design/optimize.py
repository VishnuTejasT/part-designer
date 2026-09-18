"""Multi-host, mode-driven codon optimization engine.

Implements PRODUCTION / FOLDING / BALANCED / ALL per the user's spec, with
this explicit priority when constraints conflict (never silently relaxed):

    protein identity > forbidden sites > zero rare codons (w<0.3) >
    repeats/GC/homopolymers/cryptic-motifs > CAI floor > hairpin limits

"Zero rare codons" is enforced structurally: every codon choice throughout
this module is drawn only from each amino acid's "safe" synonym set (w>=0.3
relative to that amino acid's most-used codon) -- see ``_safe_options``.
Since the most-used codon always has w=1.0, this set is never empty, so the
constraint can never actually be violated by construction, only reported.

Hairpin/secondary-structure checks use real ViennaRNA folding
(secondary_structure.py) at the caller's folding temperature, refolding
only the sliding windows that overlap a changed codon (per the spec's
"After each swap, refold only the 60 nt windows that overlap the changed
codon" instruction) rather than the whole sequence every iteration.
"""

from __future__ import annotations

import dataclasses
import random
from dataclasses import dataclass

from . import folding
from . import secondary_structure as ss
from .cai import cai_score, expression_weighted_cai_score, w_scores
from .codon_usage import CodonUsageTable, ORGANISM_TABLES, load_codon_table
from .dna_utils import gc_content, gc_window_scan, homopolymer_runs, repeated_kmers
from .forbidden_sites import scan_forbidden_sites
from .reverse_translate import (
    STOP_SYMBOL,
    CodonChoice,
    clean_protein_sequence,
    validate_protein_sequence,
)
from .sequence_motifs import scan_motifs
from .trna_usage import TrnaDataError, load_trna_table, tai_score

MODES = ("PRODUCTION", "FOLDING", "BALANCED")
DEFAULT_GC_BOUNDS = (0.30, 0.70)
RARE_CODON_W = 0.30  # hard constraint 3: no codon below this w, anywhere
CAI_FLOOR = 0.90
BALANCED_CAI_FLOOR = 0.93
PRODUCTION_CAI_TARGET = 0.97  # PRODUCTION/FOLDING maximize toward this; 0.90 stays the hard floor
PRODUCTION_DG_TARGET = 0.0  # PRODUCTION targets ~0 kcal/mol over -20/+40, not just the -5 floor
MAX_REPAIR_ITERATIONS = 300
MAX_CAI_RAISE_ITERATIONS = 150
MAX_INITIATION_ITERATIONS = 60
MAX_HAIRPIN_ITERATIONS = 60
INITIATION_WINDOW_CODONS = 13
PAUSE_W_LOW, PAUSE_W_HIGH = 0.30, 0.50
INVARIANT_AA = ("M", "W", "*")


@dataclass(frozen=True)
class StructuralRegion:
    start: int  # 1-indexed amino-acid position, inclusive
    end: int
    kind: str  # linker | loop | domain_boundary | helix | strand


@dataclass(frozen=True)
class OptimizationRequest:
    protein: str
    hosts: tuple[str, ...]
    mode: str = "BALANCED"  # PRODUCTION | FOLDING | BALANCED | ALL
    native_organism: str | None = None
    five_prime_utr: str = ""
    n_tag: str = ""
    c_tag: str = ""
    forbidden_enzymes: tuple[tuple[str, str], ...] = ()
    structural_regions: tuple[StructuralRegion, ...] = ()
    hedge: bool = False
    seed: int | None = None
    gc_bounds: tuple[float, float] = DEFAULT_GC_BOUNDS
    temperature_c: float = ss.DEFAULT_TEMPERATURE_C


@dataclass(frozen=True)
class HostModeResult:
    host: str
    mode: str
    protein: str
    dna: str
    cai: float | None
    cai_note: str | None
    expression_weighted_cai: float | None
    expression_weighted_cai_note: str | None
    tai: float | None
    tai_note: str | None
    min_w_used: float
    codons_below_w_threshold: int
    five_prime_dG: float | None
    five_prime_dG_note: str | None
    init_region_unpaired: bool | None
    worst_window_mfe: float | None
    longest_stem: int
    terminator_motifs: list
    inverted_repeats: list
    gc_overall: float
    gc_window_min: float | None
    gc_window_max: float | None
    repeated_15mers: int
    forbidden_site_hits: dict
    motif_hits: list
    pause_site_count: int
    constraint_pass_fail: dict[str, bool]
    notes: list[str]
    codon_choices: tuple[CodonChoice, ...]


class OptimizationError(ValueError):
    pass


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------


def _codon_positions_overlapping(start: int, end: int) -> list[int]:
    first = (start - 1) // 3 + 1
    last = (end - 1) // 3 + 1
    return list(range(first, last + 1))


def _safe_options(table: CodonUsageTable, amino_acid: str, scores: dict[str, float]):
    """Synonymous codons for ``amino_acid`` with w >= RARE_CODON_W. Never
    empty: the max-fraction codon always has w=1.0."""

    options = table.codons_for(amino_acid)
    safe = [o for o in options if scores.get(o.codon, 0.0) >= RARE_CODON_W]
    return safe or [max(options, key=lambda o: o.fraction)]


def _weighted_pick(options, rng: random.Random):
    weights = [o.fraction for o in options]
    if sum(weights) <= 0:
        weights = [1.0] * len(options)
    return rng.choices(options, weights=weights, k=1)[0]


def _initial_translation(
    protein: str, table: CodonUsageTable, scores: dict[str, float]
) -> tuple[list[str], list[CodonChoice]]:
    """Deterministic highest-weight codon at every position, per the spec's
    Method section. Since w=1.0 codons are used throughout, this starts at
    CAI=1.0 and trivially satisfies the zero-rare-codon constraint."""

    residues = protein if protein.endswith(STOP_SYMBOL) else protein + STOP_SYMBOL
    codons: list[str] = []
    choices: list[CodonChoice] = []
    for position, amino_acid in enumerate(residues, start=1):
        best = max(table.codons_for(amino_acid), key=lambda o: o.fraction)
        codons.append(best.codon)
        choices.append(
            CodonChoice(position=position, amino_acid=amino_acid, codon=best.codon, fraction=best.fraction)
        )
    return codons, choices


def _harmonized_translation(
    protein: str,
    host_table: CodonUsageTable,
    native_table: CodonUsageTable,
    host_scores: dict[str, float],
    rng: random.Random,
) -> tuple[list[str], list[CodonChoice]]:
    """%MinMax-style harmonization approximation: since we only have the
    protein sequence (no actual native CDS was supplied), a plausible native
    codon per position is sampled from the native table's own usage
    distribution, its usage-rank percentile within its synonymous group is
    computed, and the host codon at the matching percentile among the
    host's *safe* (w>=0.3) synonyms is chosen -- harmonization never drops
    below w=0.3, per the spec."""

    residues = protein if protein.endswith(STOP_SYMBOL) else protein + STOP_SYMBOL
    codons: list[str] = []
    choices: list[CodonChoice] = []
    for position, amino_acid in enumerate(residues, start=1):
        native_options = sorted(native_table.codons_for(amino_acid), key=lambda o: -o.fraction)
        native_pick = _weighted_pick(native_options, rng)
        native_rank = next(i for i, o in enumerate(native_options) if o.codon == native_pick.codon)
        percentile = native_rank / max(1, len(native_options) - 1)

        host_options = sorted(_safe_options(host_table, amino_acid, host_scores), key=lambda o: -o.fraction)
        host_rank = round(percentile * (len(host_options) - 1))
        chosen = host_options[host_rank]
        codons.append(chosen.codon)
        choices.append(
            CodonChoice(position=position, amino_acid=amino_acid, codon=chosen.codon, fraction=chosen.fraction)
        )
    return codons, choices


def _sequence_violations(
    dna: str,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
) -> list[tuple[int, int, str]]:
    """Priority-tier-2/3 violations: forbidden sites, homopolymers, repeats,
    GC windows, cryptic/SD-like motifs. (Rare codons are excluded here --
    they're prevented structurally, not detected-and-repaired.)"""

    violations: list[tuple[int, int, str]] = []

    for name, occs in scan_forbidden_sites(dna, forbidden_extra).items():
        for occ in occs:
            violations.append((occ.start, occ.end, f"forbidden site {name}"))

    for run in homopolymer_runs(dna):
        violations.append((run.start, run.end, "homopolymer run"))

    for rep in repeated_kmers(dna):
        # Only flag every occurrence *after* the first: if a repair pass
        # resampled every occurrence of a duplicate simultaneously, two
        # occurrences that are still identical in amino acid sequence (so
        # every candidate codon set is identical too) can flip in perfect
        # lockstep forever. Leaving the first copy untouched breaks that
        # symmetry deterministically.
        for pos in rep.positions[1:]:
            violations.append((pos, pos + 14, "repeated 15-mer"))

    low, high = gc_bounds
    for gcv in gc_window_scan(dna, low=low, high=high):
        violations.append((gcv.start, gcv.end, f"GC window {gcv.bound}"))

    if domain:
        for hit in scan_motifs(dna, domain):
            violations.append((hit.start, hit.end, f"cryptic motif ({hit.kind})"))

    return violations


def _repair_sequence_constraints(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    rng: random.Random,
    max_iterations: int = MAX_REPAIR_ITERATIONS,
) -> tuple[bool, int, list[str]]:
    """Fix forbidden-site/repeat/GC/homopolymer/motif violations, resampling
    only from each position's safe (w>=0.3) synonym set -- so a fix can
    never introduce a rare codon."""

    iterations = 0
    last_reasons: list[str] = []
    while iterations < max_iterations:
        dna = "".join(codons)
        violations = _sequence_violations(dna, domain, forbidden_extra, gc_bounds)
        if not violations:
            return True, iterations, []
        last_reasons = sorted({v[2] for v in violations})
        iterations += 1

        position_reasons: dict[int, str] = {}
        for start, end, reason in violations:
            for pos in _codon_positions_overlapping(start, end):
                if 1 <= pos <= len(codons):
                    position_reasons.setdefault(pos, reason)

        changed_any = False
        for pos in position_reasons:
            amino_acid = choices[pos - 1].amino_acid
            if amino_acid in INVARIANT_AA:
                continue
            safe = _safe_options(table, amino_acid, scores)
            current = codons[pos - 1]
            alternatives = [o for o in safe if o.codon != current]
            if not alternatives:
                continue
            chosen = _weighted_pick(alternatives, rng)
            codons[pos - 1] = chosen.codon
            choices[pos - 1] = dataclasses.replace(choices[pos - 1], codon=chosen.codon, fraction=chosen.fraction)
            changed_any = True

        if not changed_any:
            return False, iterations, last_reasons

    return False, iterations, last_reasons


def _raise_cai(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    target: float,
    max_iterations: int = MAX_CAI_RAISE_ITERATIONS,
) -> float:
    """Push CAI toward ``target`` without breaking tier-2/3 sequence
    constraints (higher priority than CAI). Only used when the initial
    max-w start got pulled down by the constraint-repair pass."""

    def current_cai() -> float:
        try:
            return cai_score("".join(codons), table)
        except ValueError:
            return 1.0

    cai = current_cai()
    # Compare against the *current* violation count, not "zero violations":
    # an unresolved tier-2 violation elsewhere in the sequence (e.g. a
    # repeat that couldn't be broken within the safe-codon set) must not
    # block every subsequent swap that doesn't itself add a new violation.
    baseline_count = len(_sequence_violations("".join(codons), domain, forbidden_extra, gc_bounds))
    iterations = 0
    while cai < target and iterations < max_iterations:
        iterations += 1
        candidates = sorted(
            ((i, scores.get(codons[i], 0.0)) for i in range(len(codons)) if choices[i].amino_acid not in INVARIANT_AA),
            key=lambda t: t[1],
        )
        improved = False
        for i, _w in candidates:
            amino_acid = choices[i].amino_acid
            best = max(_safe_options(table, amino_acid, scores), key=lambda o: o.fraction)
            if best.codon == codons[i]:
                continue
            old_codon = codons[i]
            codons[i] = best.codon
            dna = "".join(codons)
            if len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds)) > baseline_count:
                codons[i] = old_codon
                continue
            choices[i] = dataclasses.replace(choices[i], codon=best.codon, fraction=best.fraction)
            improved = True
            cai = current_cai()
            break
        if not improved:
            break

    return cai


# --------------------------------------------------------------------------
# hairpin / secondary-structure repair (lowest priority tier)
# --------------------------------------------------------------------------


@dataclass
class _HairpinState:
    windows: list  # list[ss.WindowFold]
    worst_mfe: float | None
    longest_stem: int
    terminator_hits: list
    inverted_reps: list
    init_dG: float | None
    init_unpaired: bool | None


def _scan_hairpins(
    dna: str, domain: str | None, five_prime_utr: str, temperature_c: float
) -> _HairpinState | None:
    windows = ss.sliding_window_fold(dna, temperature_c)
    if windows is None:
        return None
    worst_mfe, longest_stem = ss.worst_window_metrics(windows)
    term_hits = ss.terminator_like_motifs(dna, windows) if domain == "prokaryote" else []
    inv_reps = ss.inverted_repeats(dna)
    init_dG = folding.initiation_dG(dna, five_prime_utr, temperature_c) if five_prime_utr else None
    init_unpaired = (
        folding.initiation_unpaired_ok(dna, five_prime_utr, temperature_c) if five_prime_utr else None
    )
    return _HairpinState(
        windows=windows, worst_mfe=worst_mfe, longest_stem=longest_stem,
        terminator_hits=term_hits, inverted_reps=inv_reps, init_dG=init_dG, init_unpaired=init_unpaired,
    )


def _refold_overlapping(
    dna: str, windows: list, change_start: int, change_end: int, temperature_c: float
) -> list | None:
    """Re-fold only the windows whose span overlaps [change_start, change_end]
    (1-indexed inclusive DNA coordinates) -- everything else is untouched,
    per the spec's incremental-refold instruction."""

    updated = []
    for w in windows:
        if w.end < change_start or w.start > change_end:
            updated.append(w)
            continue
        chunk = dna[w.start - 1 : w.end]
        result = ss.fold(chunk, temperature_c)
        if result is None:
            return None
        structure, mfe_val = result
        updated.append(ss.WindowFold(start=w.start, end=w.end, structure=structure, mfe=mfe_val))
    return updated


def _hairpin_score(state: _HairpinState) -> float:
    score = 0.0
    if state.worst_mfe is not None and state.worst_mfe < ss.MAX_WINDOW_MFE:
        score += ss.MAX_WINDOW_MFE - state.worst_mfe
    if state.longest_stem >= ss.MAX_STEM_BP:
        score += state.longest_stem - ss.MAX_STEM_BP + 1
    score += len(state.terminator_hits) * 5
    score += len(state.inverted_reps) * 2
    if state.init_unpaired is False:
        score += 10
    if state.init_dG is not None and state.init_dG < ss.INIT_DG_FLOOR:
        score += ss.INIT_DG_FLOOR - state.init_dG
    return score


def _repair_hairpins(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    five_prime_utr: str,
    temperature_c: float,
    cai_floor: float,
    rng: random.Random,
    target_dG: float | None = None,
    positions: list[int] | None = None,
    max_iterations: int = MAX_HAIRPIN_ITERATIONS,
) -> tuple[_HairpinState | None, bool, bool]:
    """Lowest-priority-tier local search: try to resolve hairpin violations
    (6a-6e) via safe-synonym swaps, incrementally refolding only the
    affected windows. Never breaks tier-2/3 sequence constraints; never
    drops CAI below ``cai_floor`` (CAI floor outranks hairpin limits).

    Returns (final hairpin state, fully_resolved, blocked_by_cai_floor).
    """

    dna = "".join(codons)
    state = _scan_hairpins(dna, domain, five_prime_utr, temperature_c)
    if state is None:
        return None, False, False

    def resolved(s: _HairpinState) -> bool:
        ok = (
            (s.worst_mfe is None or s.worst_mfe >= ss.MAX_WINDOW_MFE)
            and s.longest_stem < ss.MAX_STEM_BP
            and not s.terminator_hits
            and not s.inverted_reps
            and s.init_unpaired is not False
            and (s.init_dG is None or s.init_dG >= ss.INIT_DG_FLOOR)
        )
        if target_dG is not None and s.init_dG is not None:
            ok = ok and s.init_dG >= target_dG
        return ok

    candidate_positions = positions or list(range(len(codons)))
    blocked_by_floor = False
    score = _hairpin_score(state)
    baseline_count = len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds))

    for _ in range(max_iterations):
        if resolved(state) and target_dG is None:
            break
        if not candidate_positions:
            break
        pos = rng.choice(candidate_positions)
        amino_acid = choices[pos].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        safe = _safe_options(table, amino_acid, scores)
        current = codons[pos]
        alternatives = [o for o in safe if o.codon != current]
        if not alternatives:
            continue
        candidate = rng.choice(alternatives)

        old_codon = codons[pos]
        codons[pos] = candidate.codon
        new_dna = "".join(codons)
        change_start, change_end = pos * 3 + 1, pos * 3 + 3

        if len(_sequence_violations(new_dna, domain, forbidden_extra, gc_bounds)) > baseline_count:
            codons[pos] = old_codon
            continue

        try:
            new_cai = cai_score(new_dna, table)
        except ValueError:
            new_cai = 1.0
        if new_cai < cai_floor:
            codons[pos] = old_codon
            blocked_by_floor = True
            continue

        new_windows = _refold_overlapping(new_dna, state.windows, change_start, change_end, temperature_c)
        if new_windows is None:
            codons[pos] = old_codon
            break  # ViennaRNA stopped being available mid-search
        worst_mfe, longest_stem = ss.worst_window_metrics(new_windows)
        term_hits = ss.terminator_like_motifs(new_dna, new_windows) if domain == "prokaryote" else []
        inv_reps = ss.inverted_repeats(new_dna)
        init_dG = folding.initiation_dG(new_dna, five_prime_utr, temperature_c) if five_prime_utr else None
        init_unpaired = (
            folding.initiation_unpaired_ok(new_dna, five_prime_utr, temperature_c) if five_prime_utr else None
        )
        new_state = _HairpinState(
            windows=new_windows, worst_mfe=worst_mfe, longest_stem=longest_stem,
            terminator_hits=term_hits, inverted_reps=inv_reps, init_dG=init_dG, init_unpaired=init_unpaired,
        )
        new_score = _hairpin_score(new_state)

        # For PRODUCTION's beyond-the-floor dG target, also reward pushing
        # init dG toward target_dG even once hard hairpin limits are met.
        target_bonus = 0.0
        if target_dG is not None and init_dG is not None and state.init_dG is not None:
            target_bonus = (target_dG - init_dG) - (target_dG - state.init_dG) if init_dG > state.init_dG else 0.0

        if new_score <= score and (new_score < score or target_bonus > 0):
            choices[pos] = dataclasses.replace(choices[pos], codon=candidate.codon, fraction=candidate.fraction)
            state = new_state
            score = new_score
        else:
            codons[pos] = old_codon

    return state, resolved(state), blocked_by_floor


# --------------------------------------------------------------------------
# rare-codon / w reporting (structural guarantee verified, not enforced here)
# --------------------------------------------------------------------------


def _w_stats(codons: list[str], choices: list[CodonChoice], scores: dict[str, float]) -> tuple[float, int]:
    ws = [
        scores.get(codon, 1.0) if choices[i].amino_acid not in INVARIANT_AA else 1.0
        for i, codon in enumerate(codons)
    ]
    below = sum(1 for w in ws if w < RARE_CODON_W)
    return (min(ws) if ws else 1.0), below


# --------------------------------------------------------------------------
# initiation-region search (PRODUCTION)
# --------------------------------------------------------------------------


def _optimize_initiation(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    five_prime_utr: str,
    temperature_c: float,
    floor: float,
    rng: random.Random,
    iterations: int = MAX_INITIATION_ITERATIONS,
) -> float | None:
    dna = "".join(codons)
    best_dG = folding.initiation_dG(dna, five_prime_utr, temperature_c)
    if best_dG is None:
        return None
    baseline_count = len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds))

    limit = min(INITIATION_WINDOW_CODONS, len(codons) - 1)
    for _ in range(iterations):
        if best_dG >= PRODUCTION_DG_TARGET or limit < 1:
            break
        pos = rng.randrange(1, limit + 1)
        amino_acid = choices[pos].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        safe = _safe_options(table, amino_acid, scores)
        alternatives = [o for o in safe if o.codon != codons[pos]]
        if not alternatives:
            continue
        candidate = rng.choice(alternatives)
        old_codon = codons[pos]
        codons[pos] = candidate.codon
        dna = "".join(codons)
        if len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds)) > baseline_count:
            codons[pos] = old_codon
            continue
        try:
            new_cai = cai_score(dna, table)
        except ValueError:
            new_cai = 1.0
        if new_cai < floor:
            codons[pos] = old_codon
            continue
        new_dG = folding.initiation_dG(dna, five_prime_utr, temperature_c)
        if new_dG is not None and new_dG > best_dG:
            best_dG = new_dG
            choices[pos] = dataclasses.replace(choices[pos], codon=candidate.codon, fraction=candidate.fraction)
        else:
            codons[pos] = old_codon

    return best_dG


def _insert_pause_sites(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    regions: tuple[StructuralRegion, ...],
    cai_floor: float,
    five_prime_utr: str,
    temperature_c: float,
    baseline_dG: float | None,
) -> int:
    pause_count = 0
    positions: set[int] = set()
    for region in regions:
        if region.kind not in ("linker", "loop", "domain_boundary"):
            continue
        positions.update(range(region.start, region.end + 1))

    baseline_count = len(_sequence_violations("".join(codons), domain, forbidden_extra, gc_bounds))

    for pos in sorted(positions):
        idx = pos - 1
        if idx < 0 or idx >= len(codons):
            continue
        amino_acid = choices[idx].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        options = table.codons_for(amino_acid)
        max_fraction = max(o.fraction for o in options)
        if max_fraction <= 0:
            continue
        candidates = [
            o for o in options
            if PAUSE_W_LOW <= (o.fraction / max_fraction) < PAUSE_W_HIGH
        ]
        old_codon = codons[idx]
        placed = False
        for candidate in candidates:
            codons[idx] = candidate.codon
            dna = "".join(codons)
            if len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds)) > baseline_count:
                codons[idx] = old_codon
                continue
            try:
                new_cai = cai_score(dna, table)
            except ValueError:
                new_cai = 1.0
            if new_cai < cai_floor:
                codons[idx] = old_codon
                continue
            if five_prime_utr and baseline_dG is not None:
                new_dG = folding.initiation_dG(dna, five_prime_utr, temperature_c)
                if new_dG is not None and new_dG < baseline_dG:
                    codons[idx] = old_codon
                    continue
            choices[idx] = dataclasses.replace(choices[idx], codon=candidate.codon, fraction=candidate.fraction)
            placed = True
            break
        if placed:
            pause_count += 1
        else:
            codons[idx] = old_codon

    return pause_count


def _diversify(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    other_variants: list[list[str]],
) -> None:
    baseline_count = len(_sequence_violations("".join(codons), domain, forbidden_extra, gc_bounds))
    for i in range(len(codons)):
        amino_acid = choices[i].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        safe = _safe_options(table, amino_acid, scores)
        used_elsewhere = {variant[i] for variant in other_variants}
        alternatives = [o for o in safe if o.codon not in used_elsewhere]
        if not alternatives:
            continue
        alternatives.sort(key=lambda o: -o.fraction)
        old_codon = codons[i]
        for candidate in alternatives:
            codons[i] = candidate.codon
            dna = "".join(codons)
            if len(_sequence_violations(dna, domain, forbidden_extra, gc_bounds)) > baseline_count:
                codons[i] = old_codon
                continue
            try:
                new_cai = cai_score(dna, table)
            except ValueError:
                new_cai = 1.0
            if new_cai < CAI_FLOOR:
                codons[i] = old_codon
                continue
            choices[i] = dataclasses.replace(choices[i], codon=candidate.codon, fraction=candidate.fraction)
            break
        else:
            codons[i] = old_codon


# --------------------------------------------------------------------------
# per host/mode entry point
# --------------------------------------------------------------------------


def _build_result(
    host: str,
    mode: str,
    protein: str,
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    scores: dict[str, float],
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    resolved: bool,
    blocking_reasons: list[str],
    five_prime_utr: str,
    temperature_c: float,
    hairpin_state: _HairpinState | None,
    hairpin_resolved: bool,
    hairpin_blocked_by_floor: bool,
    pause_site_count: int,
    notes: list[str],
) -> HostModeResult:
    dna = "".join(codons)

    try:
        cai = cai_score(dna, table)
        cai_note = None
    except ValueError as exc:
        cai = None
        cai_note = str(exc)

    ew_cai, ew_note = expression_weighted_cai_score(dna, table)

    tai = None
    tai_note = None
    try:
        trna_table = load_trna_table(host)
        tai = tai_score(dna, trna_table, table)
    except (TrnaDataError, ValueError) as exc:
        tai_note = f"tAI unavailable: {exc}"

    min_w, below_threshold = _w_stats(codons, choices, scores)

    gc_windows = gc_window_scan(dna, low=gc_bounds[0], high=gc_bounds[1])
    if len(dna) >= 50:
        window_fractions = [gc_content(dna[i : i + 50]) for i in range(0, len(dna) - 49)]
        gc_min, gc_max = min(window_fractions), max(window_fractions)
    else:
        gc_min = gc_max = gc_content(dna)

    forbidden_hits = scan_forbidden_sites(dna, forbidden_extra)
    motif_hits = scan_motifs(dna, domain) if domain else []
    repeats = repeated_kmers(dna)

    internal_stop_free = all(table.by_codon.get(codons[i], {}).get("amino_acid") != "*" for i in range(len(codons) - 1))
    ends_in_stop = table.by_codon.get(codons[-1], {}).get("amino_acid") == "*"

    hairpins_available = hairpin_state is not None
    constraint_pass_fail = {
        "no_internal_stop": internal_stop_free and ends_in_stop,
        "no_forbidden_sites": not forbidden_hits,
        "zero_rare_codons": below_threshold == 0,
        "no_repeated_15mers": not repeats,
        "no_long_homopolymers": not homopolymer_runs(dna),
        "gc_windows_in_bounds": not gc_windows,
        "cai_floor_met": cai is not None and cai >= CAI_FLOOR,
        "init_region_unpaired": (hairpin_state.init_unpaired is not False) if hairpins_available else True,
        "no_overlong_stem": (hairpin_state.longest_stem < ss.MAX_STEM_BP) if hairpins_available else True,
        "no_overstable_window": (
            hairpin_state.worst_mfe is None or hairpin_state.worst_mfe >= ss.MAX_WINDOW_MFE
        ) if hairpins_available else True,
        "no_terminator_motifs": (not hairpin_state.terminator_hits) if hairpins_available else True,
        "no_inverted_repeats": (not hairpin_state.inverted_reps) if hairpins_available else True,
    }

    all_notes = list(notes)
    if not resolved:
        all_notes.append(
            f"best achievable sequence: hard sequence constraints not fully resolved, "
            f"blocked by: {', '.join(blocking_reasons) or 'unknown'}"
        )
    if cai is not None and cai < CAI_FLOOR:
        all_notes.append(f"CAI floor 0.90 not met (achieved {cai:.3f})")
    if tai_note:
        all_notes.append(tai_note)
    if ew_note:
        all_notes.append(f"expression-weighted CAI unavailable: {ew_note}")

    if not hairpins_available:
        all_notes.append("hairpin/secondary-structure checks unavailable: ViennaRNA could not be loaded")
    elif not hairpin_resolved:
        parts = []
        if hairpin_state.worst_mfe is not None and hairpin_state.worst_mfe < ss.MAX_WINDOW_MFE:
            parts.append(f"worst 60nt window MFE {hairpin_state.worst_mfe:.2f} kcal/mol (limit {ss.MAX_WINDOW_MFE})")
        if hairpin_state.longest_stem >= ss.MAX_STEM_BP:
            parts.append(f"longest stem {hairpin_state.longest_stem} bp (limit <{ss.MAX_STEM_BP})")
        if hairpin_state.terminator_hits:
            parts.append(f"{len(hairpin_state.terminator_hits)} terminator-like motif(s)")
        if hairpin_state.inverted_reps:
            parts.append(f"{len(hairpin_state.inverted_reps)} inverted repeat(s) >=8bp")
        if hairpin_state.init_unpaired is False:
            parts.append("initiation region (-15/+20) not fully single-stranded")
        if hairpin_state.init_dG is not None and hairpin_state.init_dG < ss.INIT_DG_FLOOR:
            parts.append(f"initiation dG {hairpin_state.init_dG:.2f} kcal/mol (floor {ss.INIT_DG_FLOOR})")
        detail = "; ".join(parts) or "unspecified"
        if hairpin_blocked_by_floor:
            all_notes.append(
                f"hairpin limit(s) not fully met: {detail}. Blocked by the CAI floor "
                f"(0.90): further hairpin-reducing swaps would drop CAI below the floor. "
                "Choose whether to lower the CAI floor (e.g. to 0.85) or loosen the hairpin limit."
            )
        else:
            all_notes.append(
                f"hairpin limit(s) not fully met after the local-search budget: {detail}"
            )

    return HostModeResult(
        host=host, mode=mode, protein=protein, dna=dna,
        cai=cai, cai_note=cai_note,
        expression_weighted_cai=ew_cai, expression_weighted_cai_note=ew_note,
        tai=tai, tai_note=tai_note,
        min_w_used=min_w, codons_below_w_threshold=below_threshold,
        five_prime_dG=hairpin_state.init_dG if hairpins_available else None,
        five_prime_dG_note=None if (hairpins_available and hairpin_state.init_dG is not None) else (
            "5' dG unavailable: no vector 5' UTR supplied" if not five_prime_utr else
            "5' dG unavailable: ViennaRNA could not be loaded"
        ),
        init_region_unpaired=hairpin_state.init_unpaired if hairpins_available else None,
        worst_window_mfe=hairpin_state.worst_mfe if hairpins_available else None,
        longest_stem=hairpin_state.longest_stem if hairpins_available else 0,
        terminator_motifs=[dataclasses.asdict(h) for h in hairpin_state.terminator_hits] if hairpins_available else [],
        inverted_repeats=[dataclasses.asdict(h) for h in hairpin_state.inverted_reps] if hairpins_available else [],
        gc_overall=gc_content(dna), gc_window_min=gc_min, gc_window_max=gc_max,
        repeated_15mers=len(repeats),
        forbidden_site_hits={k: len(v) for k, v in forbidden_hits.items()},
        motif_hits=[dataclasses.asdict(h) for h in motif_hits],
        pause_site_count=pause_site_count,
        constraint_pass_fail=constraint_pass_fail,
        notes=all_notes,
        codon_choices=tuple(choices),
    )


def _run_production(
    protein: str, host: str, table: CodonUsageTable, scores: dict[str, float],
    domain: str | None, request: OptimizationRequest, rng: random.Random,
) -> HostModeResult:
    codons, choices = _initial_translation(protein, table, scores)
    resolved, _iters, reasons = _repair_sequence_constraints(
        codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds, rng,
    )
    _raise_cai(codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds, PRODUCTION_CAI_TARGET)

    if request.five_prime_utr:
        _optimize_initiation(
            codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds,
            request.five_prime_utr, request.temperature_c, CAI_FLOOR, rng,
        )

    limit = min(INITIATION_WINDOW_CODONS, len(codons) - 1)
    hairpin_state, hairpin_resolved, blocked_by_floor = _repair_hairpins(
        codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds,
        request.five_prime_utr, request.temperature_c, CAI_FLOOR, rng,
        target_dG=PRODUCTION_DG_TARGET if request.five_prime_utr else None,
    )

    return _build_result(
        host, "PRODUCTION", protein, codons, choices, table, scores, domain,
        request.forbidden_enzymes, request.gc_bounds, resolved, reasons,
        request.five_prime_utr, request.temperature_c,
        hairpin_state, hairpin_resolved, blocked_by_floor, 0, [],
    )


def _run_folding(
    protein: str, host: str, table: CodonUsageTable, scores: dict[str, float],
    domain: str | None, request: OptimizationRequest, rng: random.Random,
) -> HostModeResult:
    notes = []
    native_table = None
    if request.native_organism:
        if request.native_organism in ORGANISM_TABLES:
            native_table = load_codon_table(request.native_organism)
        else:
            notes.append(
                f"harmonization skipped: native_organism {request.native_organism!r} "
                "is not in the tracked chassis set, used a CAI-natural profile instead"
            )

    if native_table is not None:
        codons, choices = _harmonized_translation(protein, table, native_table, scores, rng)
        notes.append(
            "harmonization approximated from the native organism's codon-usage "
            "distribution (no native CDS was supplied, only the protein sequence)"
        )
    else:
        codons, choices = _initial_translation(protein, table, scores)
        if not request.native_organism:
            notes.append(
                "de novo protein: no native organism given, used w>=0.5-preserving "
                "natural profile instead of harmonization"
            )

    resolved, _iters, reasons = _repair_sequence_constraints(
        codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds, rng,
    )
    _raise_cai(codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds, CAI_FLOOR)

    pause_count = 0
    if request.structural_regions:
        pause_count = _insert_pause_sites(
            codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds,
            request.structural_regions, CAI_FLOOR, request.five_prime_utr, request.temperature_c, None,
        )
    else:
        notes.append(
            "pause-site placement skipped: no structural_regions supplied "
            "(protein secondary-structure prediction is not available in this pipeline)"
        )

    hairpin_state, hairpin_resolved, blocked_by_floor = _repair_hairpins(
        codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds,
        request.five_prime_utr, request.temperature_c, CAI_FLOOR, rng,
    )

    return _build_result(
        host, "FOLDING", protein, codons, choices, table, scores, domain,
        request.forbidden_enzymes, request.gc_bounds, resolved, reasons,
        request.five_prime_utr, request.temperature_c,
        hairpin_state, hairpin_resolved, blocked_by_floor, pause_count, notes,
    )


def _run_balanced(
    protein: str, host: str, table: CodonUsageTable, scores: dict[str, float],
    domain: str | None, request: OptimizationRequest, rng: random.Random,
) -> HostModeResult:
    production = _run_production(protein, host, table, scores, domain, request, rng)
    codons = [production.dna[i : i + 3] for i in range(0, len(production.dna), 3)]
    choices = list(production.codon_choices)
    baseline_dG = production.five_prime_dG

    notes = list(production.notes)
    pause_count = 0
    if request.structural_regions:
        pause_count = _insert_pause_sites(
            codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds,
            request.structural_regions, BALANCED_CAI_FLOOR, request.five_prime_utr, request.temperature_c, baseline_dG,
        )
    else:
        notes.append(
            "pause-site placement skipped: no structural_regions supplied "
            "(protein secondary-structure prediction is not available in this pipeline)"
        )

    dna = "".join(codons)
    hairpin_state = _scan_hairpins(dna, domain, request.five_prime_utr, request.temperature_c)
    # BALANCED doesn't re-run the hairpin search (PRODUCTION already did, and
    # pause-site insertion above is already gated on not worsening init dG);
    # just re-scan to report the final state honestly.
    if hairpin_state is not None:
        resolved_flag = (
            (hairpin_state.worst_mfe is None or hairpin_state.worst_mfe >= ss.MAX_WINDOW_MFE)
            and hairpin_state.longest_stem < ss.MAX_STEM_BP
            and not hairpin_state.terminator_hits
            and not hairpin_state.inverted_reps
            and hairpin_state.init_unpaired is not False
            and (hairpin_state.init_dG is None or hairpin_state.init_dG >= ss.INIT_DG_FLOOR)
        )
    else:
        resolved_flag = False

    return _build_result(
        host, "BALANCED", protein, codons, choices, table, scores, domain,
        request.forbidden_enzymes, request.gc_bounds, True, [],
        request.five_prime_utr, request.temperature_c,
        hairpin_state, resolved_flag, False, pause_count, notes,
    )


def optimize_cds(request: OptimizationRequest) -> list[HostModeResult]:
    protein = request.n_tag + clean_protein_sequence(request.protein) + request.c_tag
    validate_protein_sequence(protein)

    unknown_hosts = [h for h in request.hosts if h not in ORGANISM_TABLES]
    if unknown_hosts:
        available = ", ".join(sorted(ORGANISM_TABLES))
        raise OptimizationError(f"Unknown host(s) {unknown_hosts}. Available: {available}")

    modes = MODES if request.mode == "ALL" else (request.mode,)
    if any(m not in MODES for m in modes):
        raise OptimizationError(f"Unknown mode {request.mode!r}. Choose from ALL, {', '.join(MODES)}")

    results: list[HostModeResult] = []
    for host in request.hosts:
        table = load_codon_table(host)
        scores = w_scores(table)
        domain = table.metadata.get("domain")
        base_seed = request.seed
        host_variants: list[list[str]] = []
        balanced_result = None

        for mode in modes:
            seed = None if base_seed is None else base_seed + hash((host, mode)) % 10_000
            rng = random.Random(seed)
            if mode == "PRODUCTION":
                result = _run_production(protein, host, table, scores, domain, request, rng)
            elif mode == "FOLDING":
                result = _run_folding(protein, host, table, scores, domain, request, rng)
            else:
                result = _run_balanced(protein, host, table, scores, domain, request, rng)
                balanced_result = result
            results.append(result)
            host_variants.append([result.dna[i : i + 3] for i in range(0, len(result.dna), 3)])

        if request.hedge:
            if balanced_result is None:
                seed = None if base_seed is None else base_seed + hash((host, "BALANCED")) % 10_000
                rng = random.Random(seed)
                balanced_result = _run_balanced(protein, host, table, scores, domain, request, rng)
            codons = [balanced_result.dna[i : i + 3] for i in range(0, len(balanced_result.dna), 3)]
            choices = list(balanced_result.codon_choices)
            _diversify(codons, choices, table, scores, domain, request.forbidden_enzymes, request.gc_bounds, host_variants)

            dna = "".join(codons)
            hairpin_state = _scan_hairpins(dna, domain, request.five_prime_utr, request.temperature_c)
            if hairpin_state is not None:
                resolved_flag = (
                    (hairpin_state.worst_mfe is None or hairpin_state.worst_mfe >= ss.MAX_WINDOW_MFE)
                    and hairpin_state.longest_stem < ss.MAX_STEM_BP
                    and not hairpin_state.terminator_hits
                    and not hairpin_state.inverted_reps
                    and hairpin_state.init_unpaired is not False
                    and (hairpin_state.init_dG is None or hairpin_state.init_dG >= ss.INIT_DG_FLOOR)
                )
            else:
                resolved_flag = False

            hedge_result = _build_result(
                host, "BALANCED (diversified)", protein, codons, choices, table, scores, domain,
                request.forbidden_enzymes, request.gc_bounds, True, [],
                request.five_prime_utr, request.temperature_c,
                hairpin_state, resolved_flag, False, balanced_result.pause_site_count,
                list(balanced_result.notes) + ["HEDGE variant: diversified from other mode results at CAI>=0.90"],
            )
            results.append(hedge_result)

    return results
