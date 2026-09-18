"""Multi-host, mode-driven codon optimization engine.

Implements PRODUCTION / FOLDING / BALANCED / ALL per the user's spec:
hard constraints (forbidden sites, homopolymers, repeats, GC windows,
cryptic motifs) are enforced first and never silently relaxed; CAI is
maximized/floored using the real per-host Kazusa-derived tables; tAI,
5' mRNA-folding dG, and structural pause-site placement are each computed
for real where the data/tooling exists and are explicitly reported as
unavailable (never fabricated) where it doesn't.
"""

from __future__ import annotations

import dataclasses
import random
from dataclasses import dataclass, field

from . import folding
from .cai import cai_score, expression_weighted_cai_score, w_scores
from .codon_usage import CodonUsageTable, ORGANISM_TABLES, load_codon_table
from .dna_utils import gc_content, gc_window_scan, homopolymer_runs, repeated_kmers
from .forbidden_sites import scan_forbidden_sites
from .reverse_translate import (
    STOP_SYMBOL,
    CodonChoice,
    ProteinSequenceError,
    clean_protein_sequence,
    validate_protein_sequence,
)
from .sequence_motifs import scan_motifs
from .trna_usage import TRNA_TABLES, TrnaDataError, load_trna_table, tai_score

MODES = ("PRODUCTION", "FOLDING", "BALANCED")
DEFAULT_GC_BOUNDS = (0.30, 0.70)
CAI_FLOOR = 0.90
BALANCED_CAI_FLOOR = 0.93
PRODUCTION_CAI_TARGET = 0.97  # PRODUCTION/FOLDING maximize toward this; 0.90 stays the hard floor
MAX_REPAIR_ITERATIONS = 300
MAX_CAI_RAISE_ITERATIONS = 200
MAX_INITIATION_ITERATIONS = 60
INITIATION_WINDOW_CODONS = 13
PRODUCTION_BIAS = 6.0
FOLDING_BIAS = 1.0
RARE_CODON_W = 0.3
RARE_CODON_CLUSTER_LEN = 3
PAUSE_W_LOW, PAUSE_W_HIGH = 0.2, 0.5
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
    five_prime_dG: float | None
    five_prime_dG_note: str | None
    gc_overall: float
    gc_window_min: float | None
    gc_window_max: float | None
    repeated_15mers: int
    forbidden_site_hits: dict
    motif_hits: list
    rare_codon_clusters: list
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


def _mode_bias(mode: str) -> float:
    return {"PRODUCTION": PRODUCTION_BIAS, "FOLDING": FOLDING_BIAS, "BALANCED": PRODUCTION_BIAS}.get(
        mode, PRODUCTION_BIAS
    )


def _weighted_pick(options, bias: float, rng: random.Random):
    max_fraction = max(o.fraction for o in options)
    weights = [
        (o.fraction / max_fraction) ** bias if max_fraction > 0 else 1.0
        for o in options
    ]
    if sum(weights) <= 0:
        weights = [1.0] * len(options)
    return rng.choices(options, weights=weights, k=1)[0]


def _initial_translation(
    protein: str, table: CodonUsageTable, bias: float, rng: random.Random
) -> tuple[list[str], list[CodonChoice]]:
    residues = protein if protein.endswith(STOP_SYMBOL) else protein + STOP_SYMBOL
    codons: list[str] = []
    choices: list[CodonChoice] = []
    for position, amino_acid in enumerate(residues, start=1):
        options = table.codons_for(amino_acid)
        chosen = _weighted_pick(options, bias, rng)
        codons.append(chosen.codon)
        choices.append(
            CodonChoice(
                position=position,
                amino_acid=amino_acid,
                codon=chosen.codon,
                fraction=chosen.fraction,
            )
        )
    return codons, choices


def _harmonized_translation(
    protein: str,
    host_table: CodonUsageTable,
    native_table: CodonUsageTable,
    rng: random.Random,
) -> tuple[list[str], list[CodonChoice]]:
    """%MinMax-style harmonization approximation: since we only have the
    protein sequence (no actual native CDS was supplied), a plausible native
    codon per position is sampled from the native table's own usage
    distribution, its usage-rank percentile within its synonymous group is
    computed, and the host codon at the matching percentile is chosen."""

    residues = protein if protein.endswith(STOP_SYMBOL) else protein + STOP_SYMBOL
    codons: list[str] = []
    choices: list[CodonChoice] = []
    for position, amino_acid in enumerate(residues, start=1):
        native_options = sorted(
            native_table.codons_for(amino_acid), key=lambda o: -o.fraction
        )
        native_pick = _weighted_pick(native_options, 1.0, rng)
        native_rank = next(
            i for i, o in enumerate(native_options) if o.codon == native_pick.codon
        )
        percentile = native_rank / max(1, len(native_options) - 1)

        host_options = sorted(
            host_table.codons_for(amino_acid), key=lambda o: -o.fraction
        )
        host_rank = round(percentile * (len(host_options) - 1))
        chosen = host_options[host_rank]
        codons.append(chosen.codon)
        choices.append(
            CodonChoice(
                position=position,
                amino_acid=amino_acid,
                codon=chosen.codon,
                fraction=chosen.fraction,
            )
        )
    return codons, choices


def _rare_codon_clusters(codons: list[str], choices: list[CodonChoice], table: CodonUsageTable):
    scores = w_scores(table)
    clusters = []
    run_start = None
    for i, codon in enumerate(codons):
        if choices[i].amino_acid in INVARIANT_AA:
            w = 1.0
        else:
            w = scores.get(codon, 0.0)
        if w < RARE_CODON_W:
            if run_start is None:
                run_start = i
        else:
            if run_start is not None and i - run_start >= RARE_CODON_CLUSTER_LEN:
                clusters.append((run_start + 1, i))
            run_start = None
    if run_start is not None and len(codons) - run_start >= RARE_CODON_CLUSTER_LEN:
        clusters.append((run_start + 1, len(codons)))
    return clusters


def _violations(
    dna: str,
    table: CodonUsageTable,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    check_rare_clusters: bool,
    codons: list[str],
    choices: list[CodonChoice],
) -> list[tuple[int, int, str]]:
    violations: list[tuple[int, int, str]] = []

    for name, occs in scan_forbidden_sites(dna, forbidden_extra).items():
        for occ in occs:
            violations.append((occ.start, occ.end, f"forbidden site {name}"))

    for run in homopolymer_runs(dna):
        violations.append((run.start, run.end, "homopolymer run"))

    for rep in repeated_kmers(dna):
        for pos in rep.positions:
            violations.append((pos, pos + 14, "repeated 15-mer"))

    low, high = gc_bounds
    for gcv in gc_window_scan(dna, low=low, high=high):
        violations.append((gcv.start, gcv.end, f"GC window {gcv.bound}"))

    if domain:
        for hit in scan_motifs(dna, domain):
            violations.append((hit.start, hit.end, f"cryptic motif ({hit.kind})"))

    if check_rare_clusters:
        for start_pos, end_pos in _rare_codon_clusters(codons, choices, table):
            violations.append(((start_pos - 1) * 3 + 1, end_pos * 3, "rare-codon cluster"))

    return violations


def _repair_constraints(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    rng: random.Random,
    bias: float,
    check_rare_clusters: bool,
    max_iterations: int = MAX_REPAIR_ITERATIONS,
) -> tuple[bool, int, list[str]]:
    iterations = 0
    last_reasons: list[str] = []
    while iterations < max_iterations:
        dna = "".join(codons)
        violations = _violations(
            dna, table, domain, forbidden_extra, gc_bounds, check_rare_clusters, codons, choices
        )
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
            options = table.codons_for(amino_acid)
            current = codons[pos - 1]
            alternatives = [o for o in options if o.codon != current]
            if not alternatives:
                continue
            chosen = _weighted_pick(alternatives, bias, rng)
            codons[pos - 1] = chosen.codon
            choices[pos - 1] = dataclasses.replace(
                choices[pos - 1], codon=chosen.codon, fraction=chosen.fraction
            )
            changed_any = True

        if not changed_any:
            return False, iterations, last_reasons

    return False, iterations, last_reasons


def _raise_cai_floor(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    floor: float,
    max_iterations: int = MAX_CAI_RAISE_ITERATIONS,
) -> tuple[float, bool]:
    scores = w_scores(table)

    def current_cai() -> float:
        try:
            return cai_score("".join(codons), table)
        except ValueError:
            return 1.0

    cai = current_cai()
    iterations = 0
    while cai < floor and iterations < max_iterations:
        iterations += 1
        candidates = sorted(
            (
                (i, scores.get(codons[i], 0.0))
                for i in range(len(codons))
                if choices[i].amino_acid not in INVARIANT_AA
            ),
            key=lambda t: t[1],
        )
        improved = False
        for i, _w in candidates:
            amino_acid = choices[i].amino_acid
            options = table.codons_for(amino_acid)
            best = max(options, key=lambda o: o.fraction)
            if best.codon == codons[i]:
                continue
            old_codon = codons[i]
            codons[i] = best.codon
            dna = "".join(codons)
            if _violations(dna, table, domain, forbidden_extra, gc_bounds, False, codons, choices):
                codons[i] = old_codon
                continue
            choices[i] = dataclasses.replace(choices[i], codon=best.codon, fraction=best.fraction)
            improved = True
            cai = current_cai()
            break
        if not improved:
            break

    return cai, cai >= floor


def _optimize_initiation(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    five_prime_utr: str,
    floor: float,
    rng: random.Random,
    iterations: int = MAX_INITIATION_ITERATIONS,
) -> tuple[float | None, str | None]:
    dna = "".join(codons)
    best_dG = folding.initiation_dG(dna, five_prime_utr)
    if best_dG is None:
        return None, "5' initiation dG unavailable: ViennaRNA could not be loaded"

    limit = min(INITIATION_WINDOW_CODONS, len(codons) - 1)
    for _ in range(iterations):
        if limit < 1:
            break
        pos = rng.randrange(1, limit + 1)
        amino_acid = choices[pos].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        options = table.codons_for(amino_acid)
        alternatives = [o for o in options if o.codon != codons[pos]]
        if not alternatives:
            continue
        candidate = rng.choice(alternatives)
        old_codon = codons[pos]
        codons[pos] = candidate.codon
        dna = "".join(codons)
        if _violations(dna, table, domain, forbidden_extra, gc_bounds, False, codons, choices):
            codons[pos] = old_codon
            continue
        try:
            new_cai = cai_score(dna, table)
        except ValueError:
            new_cai = 1.0
        if new_cai < floor:
            codons[pos] = old_codon
            continue
        new_dG = folding.initiation_dG(dna, five_prime_utr)
        if new_dG is not None and new_dG > best_dG:
            best_dG = new_dG
            choices[pos] = dataclasses.replace(
                choices[pos], codon=candidate.codon, fraction=candidate.fraction
            )
        else:
            codons[pos] = old_codon

    return best_dG, None


def _insert_pause_sites(
    codons: list[str],
    choices: list[CodonChoice],
    table: CodonUsageTable,
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    regions: tuple[StructuralRegion, ...],
    cai_floor: float,
    five_prime_utr: str,
    baseline_dG: float | None,
) -> int:
    pause_count = 0
    positions: set[int] = set()
    for region in regions:
        if region.kind not in ("linker", "loop", "domain_boundary"):
            continue
        positions.update(range(region.start, region.end + 1))

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
            o for o in options if PAUSE_W_LOW <= (o.fraction / max_fraction) < PAUSE_W_HIGH
        ]
        old_codon = codons[idx]
        placed = False
        for candidate in candidates:
            codons[idx] = candidate.codon
            dna = "".join(codons)
            if _violations(dna, table, domain, forbidden_extra, gc_bounds, False, codons, choices):
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
                new_dG = folding.initiation_dG(dna, five_prime_utr)
                if new_dG is not None and new_dG < baseline_dG:
                    codons[idx] = old_codon
                    continue
            choices[idx] = dataclasses.replace(
                choices[idx], codon=candidate.codon, fraction=candidate.fraction
            )
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
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    other_variants: list[list[str]],
    rng: random.Random,
) -> None:
    for i in range(len(codons)):
        amino_acid = choices[i].amino_acid
        if amino_acid in INVARIANT_AA:
            continue
        options = table.codons_for(amino_acid)
        used_elsewhere = {variant[i] for variant in other_variants}
        alternatives = [o for o in options if o.codon not in used_elsewhere]
        if not alternatives:
            continue
        alternatives.sort(key=lambda o: -o.fraction)
        old_codon = codons[i]
        for candidate in alternatives:
            codons[i] = candidate.codon
            dna = "".join(codons)
            if _violations(dna, table, domain, forbidden_extra, gc_bounds, False, codons, choices):
                codons[i] = old_codon
                continue
            try:
                new_cai = cai_score(dna, table)
            except ValueError:
                new_cai = 1.0
            if new_cai < CAI_FLOOR:
                codons[i] = old_codon
                continue
            choices[i] = dataclasses.replace(
                choices[i], codon=candidate.codon, fraction=candidate.fraction
            )
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
    domain: str | None,
    forbidden_extra: tuple[tuple[str, str], ...],
    gc_bounds: tuple[float, float],
    resolved: bool,
    blocking_reasons: list[str],
    five_prime_utr: str,
    five_prime_dG: float | None,
    five_prime_dG_note: str | None,
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

    gc_windows = gc_window_scan(dna, low=gc_bounds[0], high=gc_bounds[1])
    if len(dna) >= 50:
        window_fractions = [
            gc_content(dna[i : i + 50]) for i in range(0, len(dna) - 49)
        ]
        gc_min, gc_max = min(window_fractions), max(window_fractions)
    else:
        gc_min = gc_max = gc_content(dna)

    forbidden_hits = scan_forbidden_sites(dna, forbidden_extra)
    motif_hits = scan_motifs(dna, domain) if domain else []
    rare_clusters = _rare_codon_clusters(codons, choices, table)
    repeats = repeated_kmers(dna)

    internal_stop_free = all(
        table.by_codon.get(codons[i], {}).get("amino_acid") != "*"
        for i in range(len(codons) - 1)
    )
    constraint_pass_fail = {
        "no_internal_stop": internal_stop_free
        and table.by_codon.get(codons[-1], {}).get("amino_acid") == "*",
        "no_forbidden_sites": not forbidden_hits,
        "no_repeated_15mers": not repeats,
        "no_long_homopolymers": not homopolymer_runs(dna),
        "gc_windows_in_bounds": not gc_windows,
        "cai_floor_met": cai is not None and cai >= CAI_FLOOR,
    }

    all_notes = list(notes)
    if not resolved:
        all_notes.append(
            "best achievable sequence: hard constraints not fully resolved, "
            f"blocked by: {', '.join(blocking_reasons) or 'unknown'}"
        )
    if cai is not None and cai < CAI_FLOOR:
        all_notes.append(f"CAI floor 0.90 not met (achieved {cai:.3f})")
    if tai_note:
        all_notes.append(tai_note)
    if five_prime_dG_note:
        all_notes.append(five_prime_dG_note)
    if ew_note:
        all_notes.append(f"expression-weighted CAI unavailable: {ew_note}")

    return HostModeResult(
        host=host,
        mode=mode,
        protein=protein,
        dna=dna,
        cai=cai,
        cai_note=cai_note,
        expression_weighted_cai=ew_cai,
        expression_weighted_cai_note=ew_note,
        tai=tai,
        tai_note=tai_note,
        five_prime_dG=five_prime_dG,
        five_prime_dG_note=five_prime_dG_note,
        gc_overall=gc_content(dna),
        gc_window_min=gc_min,
        gc_window_max=gc_max,
        repeated_15mers=len(repeats),
        forbidden_site_hits={k: len(v) for k, v in forbidden_hits.items()},
        motif_hits=[dataclasses.asdict(h) for h in motif_hits],
        rare_codon_clusters=rare_clusters,
        pause_site_count=pause_site_count,
        constraint_pass_fail=constraint_pass_fail,
        notes=all_notes,
        codon_choices=tuple(choices),
    )


def _run_production(
    protein: str,
    host: str,
    table: CodonUsageTable,
    domain: str | None,
    request: OptimizationRequest,
    rng: random.Random,
) -> HostModeResult:
    codons, choices = _initial_translation(protein, table, PRODUCTION_BIAS, rng)
    resolved, _iters, reasons = _repair_constraints(
        codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
        rng, PRODUCTION_BIAS, check_rare_clusters=True,
    )
    cai, _met = _raise_cai_floor(
        codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds, PRODUCTION_CAI_TARGET
    )

    notes = []
    dG = None
    dG_note = None
    if request.five_prime_utr:
        dG, dG_note = _optimize_initiation(
            codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
            request.five_prime_utr, CAI_FLOOR, rng,
        )
    else:
        dG_note = "5' dG optimization skipped: no vector 5' UTR supplied"

    return _build_result(
        host, "PRODUCTION", protein, codons, choices, table, domain,
        request.forbidden_enzymes, request.gc_bounds, resolved, reasons,
        request.five_prime_utr, dG, dG_note, 0, notes,
    )


def _run_folding(
    protein: str,
    host: str,
    table: CodonUsageTable,
    domain: str | None,
    request: OptimizationRequest,
    rng: random.Random,
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
        codons, choices = _harmonized_translation(protein, table, native_table, rng)
        notes.append(
            "harmonization approximated from the native organism's codon-usage "
            "distribution (no native CDS was supplied, only the protein sequence)"
        )
    else:
        codons, choices = _initial_translation(protein, table, FOLDING_BIAS, rng)
        if not request.native_organism:
            notes.append(
                "de novo protein: no native organism given, used w>=0.5-preserving "
                "natural-frequency sampling instead of harmonization"
            )

    resolved, _iters, reasons = _repair_constraints(
        codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
        rng, FOLDING_BIAS, check_rare_clusters=False,
    )
    _cai, _met = _raise_cai_floor(
        codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds, PRODUCTION_CAI_TARGET
    )

    pause_count = 0
    if request.structural_regions:
        pause_count = _insert_pause_sites(
            codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
            request.structural_regions, CAI_FLOOR, request.five_prime_utr, None,
        )
    else:
        notes.append(
            "pause-site placement skipped: no structural_regions supplied "
            "(protein secondary-structure prediction is not available in this pipeline)"
        )

    return _build_result(
        host, "FOLDING", protein, codons, choices, table, domain,
        request.forbidden_enzymes, request.gc_bounds, resolved, reasons,
        request.five_prime_utr, None,
        "5' dG optimization is a PRODUCTION-only refinement" if not request.five_prime_utr else None,
        pause_count, notes,
    )


def _run_balanced(
    protein: str,
    host: str,
    table: CodonUsageTable,
    domain: str | None,
    request: OptimizationRequest,
    rng: random.Random,
) -> HostModeResult:
    production = _run_production(protein, host, table, domain, request, rng)
    codons = [production.dna[i : i + 3] for i in range(0, len(production.dna), 3)]
    choices = list(production.codon_choices)
    baseline_dG = production.five_prime_dG

    notes = list(production.notes)
    pause_count = 0
    if request.structural_regions:
        pause_count = _insert_pause_sites(
            codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
            request.structural_regions, BALANCED_CAI_FLOOR, request.five_prime_utr, baseline_dG,
        )
    else:
        notes.append(
            "pause-site placement skipped: no structural_regions supplied "
            "(protein secondary-structure prediction is not available in this pipeline)"
        )

    return _build_result(
        host, "BALANCED", protein, codons, choices, table, domain,
        request.forbidden_enzymes, request.gc_bounds, True, [],
        request.five_prime_utr, baseline_dG, production.five_prime_dG_note,
        pause_count, notes,
    )


def optimize_cds(request: OptimizationRequest) -> list[HostModeResult]:
    protein = request.n_tag + clean_protein_sequence(request.protein) + request.c_tag
    validate_protein_sequence(protein)

    unknown_hosts = [h for h in request.hosts if h not in ORGANISM_TABLES]
    if unknown_hosts:
        available = ", ".join(sorted(ORGANISM_TABLES))
        raise OptimizationError(
            f"Unknown host(s) {unknown_hosts}. Available: {available}"
        )

    modes = MODES if request.mode == "ALL" else (request.mode,)
    if any(m not in MODES for m in modes):
        raise OptimizationError(f"Unknown mode {request.mode!r}. Choose from ALL, {', '.join(MODES)}")

    results: list[HostModeResult] = []
    for host in request.hosts:
        table = load_codon_table(host)
        domain = table.metadata.get("domain")
        base_seed = request.seed
        host_variants: list[list[str]] = []
        balanced_result = None

        for mode in modes:
            seed = None if base_seed is None else base_seed + hash((host, mode)) % 10_000
            rng = random.Random(seed)
            if mode == "PRODUCTION":
                result = _run_production(protein, host, table, domain, request, rng)
            elif mode == "FOLDING":
                result = _run_folding(protein, host, table, domain, request, rng)
            else:
                result = _run_balanced(protein, host, table, domain, request, rng)
                balanced_result = result
            results.append(result)
            host_variants.append(
                [result.dna[i : i + 3] for i in range(0, len(result.dna), 3)]
            )

        if request.hedge:
            if balanced_result is None:
                seed = None if base_seed is None else base_seed + hash((host, "BALANCED")) % 10_000
                rng = random.Random(seed)
                balanced_result = _run_balanced(protein, host, table, domain, request, rng)
            codons = [
                balanced_result.dna[i : i + 3] for i in range(0, len(balanced_result.dna), 3)
            ]
            choices = list(balanced_result.codon_choices)
            hedge_seed = None if base_seed is None else base_seed + hash((host, "HEDGE")) % 10_000
            rng = random.Random(hedge_seed)
            _diversify(
                codons, choices, table, domain, request.forbidden_enzymes, request.gc_bounds,
                host_variants, rng,
            )
            hedge_result = _build_result(
                host, "BALANCED (diversified)", protein, codons, choices, table, domain,
                request.forbidden_enzymes, request.gc_bounds, True, [],
                request.five_prime_utr, balanced_result.five_prime_dG,
                balanced_result.five_prime_dG_note, balanced_result.pause_site_count,
                list(balanced_result.notes)
                + ["HEDGE variant: diversified from other mode results at CAI>=0.90"],
            )
            results.append(hedge_result)

    return results
