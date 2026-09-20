"""Engine fields added for the UI redesign. Defaults must not change output."""
import json
from pathlib import Path

import pytest

from plasmid_design.optimize import (
    MAX_PROTEIN_LENGTH,
    Limits,
    OptimizationError,
    OptimizationRequest,
    optimize_cds,
)
from plasmid_design.optimize_report import build_optimization_report

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
BASELINE = json.loads((Path(__file__).parent / "fixtures" / "parity_baseline.json").read_text())


@pytest.mark.parametrize("case", BASELINE, ids=lambda c: f"{c['host']}-{c['mode']}-{c['seed']}")
def test_defaults_reproduce_baseline_dna(case):
    """Acceptance #9: same request and seed => same DNA. Explicit default
    limits must be indistinguishable from not passing any."""
    req = OptimizationRequest(
        protein=case["protein"], hosts=(case["host"],), mode=case["mode"],
        five_prime_utr=case["utr"], seed=case["seed"], limits=Limits(),
    )
    assert optimize_cds(req)[0].dna == case["dna"]


def test_same_seed_is_repeatable_across_calls():
    req = OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="ALL", seed=42)
    assert [r.dna for r in optimize_cds(req)] == [r.dna for r in optimize_cds(req)]


def test_lower_cai_floor_is_applied_and_reported():
    req = OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1,
                              limits=Limits(cai_floor=0.85))
    report = build_optimization_report(optimize_cds(req), req.limits)
    assert report["limits_used"]["cai_floor"] == 0.85
    codon_score = next(c for c in report["sequences"][0]["checks"] if c["id"] == "codon_score")
    assert codon_score["value"]["goal"] == 0.85


def test_limits_do_not_leak_between_calls():
    default = OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1)
    before = optimize_cds(default)[0].dna
    optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1,
                                     limits=Limits(rare_codon_w=0.5, cai_floor=0.5)))
    assert optimize_cds(default)[0].dna == before


def test_stricter_rare_codon_cutoff_is_enforced():
    r = optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1,
                                         limits=Limits(rare_codon_w=0.5)))[0]
    assert r.min_w_used >= 0.5 and r.codons_below_w_threshold == 0


def test_vector_provides_start_drops_the_atg():
    base = optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1))[0]
    r = optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1,
                                         vector_provides_start=True))[0]
    assert base.dna.startswith("ATG") and r.dna == base.dna[3:]
    checks = build_optimization_report([r])["sequences"][0]["checks"]
    identity = next(c for c in checks if c["id"] == "protein_identity")
    assert identity["status"] == "pass"
    assert identity["value"]["matched"] == identity["value"]["total"] == len(UBIQUITIN) - 1


def test_vector_provides_start_needs_leading_met():
    with pytest.raises(OptimizationError, match="starts with M"):
        optimize_cds(OptimizationRequest(protein="KTAYIAK", hosts=("e_coli_k12",), vector_provides_start=True))


def test_length_cap_message():
    with pytest.raises(OptimizationError, match=r"The limit is 1,500\. Try splitting it into domains\."):
        optimize_cds(OptimizationRequest(protein="M" + "A" * MAX_PROTEIN_LENGTH, hosts=("e_coli_k12",)))


def test_checks_shape_and_review_when_no_start_dna():
    report = build_optimization_report(
        optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12", "human"), mode="PRODUCTION", seed=1))
    )
    ecoli, human = report["sequences"]
    ids = [c["id"] for c in ecoli["checks"]]
    assert ids == ["protein_identity", "start_stop", "cut_sites", "rare_codons", "codon_score",
                   "synthesis", "start_region", "hairpins", "terminators", "mirror_repeats"]
    assert next(c for c in ecoli["checks"] if c["id"] == "start_region")["status"] == "review"
    assert ecoli["checks_total"] == 10
    assert human["checks_total"] == 9  # terminator check is bacteria-only


def test_conflict_is_structured_when_hairpin_limit_blocked_by_cai_floor():
    req = OptimizationRequest(protein=UBIQUITIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1,
                              five_prime_utr="AGGAGGACAGCTATG")
    conflicts = build_optimization_report(optimize_cds(req))["sequences"][0]["conflicts"]
    for c in conflicts:
        assert {"constraint", "blocked_by"} <= set(c)
        assert c["blocked_by"] in {"cai_floor", "search_budget", "sequence_rules", "protein_sequence"}


def test_forbidden_site_never_survives_when_fixable():
    """Regression: an unfixable heuristic motif must not re-break a fixed cut site."""
    r = optimize_cds(OptimizationRequest(protein=UBIQUITIN, hosts=("human",), mode="BALANCED", seed=1))[0]
    assert r.forbidden_site_hits == {}
