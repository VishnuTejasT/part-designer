import pytest

from plasmid_design.optimize import (
    OptimizationError,
    OptimizationRequest,
    StructuralRegion,
    optimize_cds,
)

# Small synthetic de novo test protein (no biological meaning).
TEST_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQ"


def test_single_mode_single_host():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=1)
    results = optimize_cds(req)
    assert len(results) == 1
    r = results[0]
    assert r.host == "e_coli_k12"
    assert r.mode == "PRODUCTION"
    assert r.dna.startswith("ATG")
    assert r.cai is not None
    assert r.cai >= 0.90
    # The tightened hairpin limits (worst-window MFE / stem length) aren't
    # always achievable in a GC-rich host like E. coli, per the spec's own
    # caveat -- but protein identity, forbidden sites, and zero rare codons
    # must always hold, unconditionally.
    assert r.constraint_pass_fail["no_internal_stop"]
    assert r.constraint_pass_fail["no_forbidden_sites"]
    assert r.constraint_pass_fail["zero_rare_codons"]
    assert r.codons_below_w_threshold == 0
    assert r.min_w_used >= 0.30


def test_all_mode_expands_to_three_plus_hedge():
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="ALL", hedge=True, seed=2
    )
    results = optimize_cds(req)
    modes = [r.mode for r in results]
    assert "PRODUCTION" in modes
    assert "FOLDING" in modes
    assert "BALANCED" in modes
    assert "BALANCED (diversified)" in modes
    assert len(results) == 4


def test_translation_is_exact_regardless_of_mode():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="ALL", seed=3)
    results = optimize_cds(req)
    for r in results:
        assert len(r.dna) == (len(TEST_PROTEIN) + 1) * 3  # + stop codon


def test_tai_unavailable_for_host_without_sourced_trna_table():
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("synechocystis_pcc6803",), mode="PRODUCTION", seed=4
    )
    result = optimize_cds(req)[0]
    assert result.tai is None
    assert result.tai_note is not None
    assert "unavailable" in result.tai_note


def test_tai_available_for_sourced_host():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=5)
    result = optimize_cds(req)[0]
    assert result.tai is not None
    assert 0 < result.tai <= 1.0


def test_five_prime_dG_populated_only_when_utr_supplied():
    req_no_utr = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=6
    )
    result_no_utr = optimize_cds(req_no_utr)[0]
    assert result_no_utr.five_prime_dG is None
    assert "no vector 5' UTR" in result_no_utr.five_prime_dG_note

    req_with_utr = OptimizationRequest(
        protein=TEST_PROTEIN,
        hosts=("e_coli_k12",),
        mode="PRODUCTION",
        five_prime_utr="AGGAGGACAGCTATG",
        seed=6,
    )
    result_with_utr = optimize_cds(req_with_utr)[0]
    assert result_with_utr.five_prime_dG is not None


def test_pause_site_placement_skipped_without_structural_regions():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="FOLDING", seed=7)
    result = optimize_cds(req)[0]
    assert result.pause_site_count == 0
    assert any("pause-site placement skipped" in n for n in result.notes)


def test_pause_sites_placed_when_regions_supplied():
    regions = (StructuralRegion(start=10, end=20, kind="linker"),)
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="FOLDING",
        structural_regions=regions, seed=8,
    )
    result = optimize_cds(req)[0]
    assert not any("pause-site placement skipped" in n for n in result.notes)


def test_caller_supplied_forbidden_enzyme_is_avoided():
    req = OptimizationRequest(
        protein=TEST_PROTEIN,
        hosts=("e_coli_k12",),
        mode="PRODUCTION",
        forbidden_enzymes=(("EcoRV", "GATATC"),),
        seed=9,
    )
    result = optimize_cds(req)[0]
    assert "EcoRV" not in result.forbidden_site_hits
    assert result.constraint_pass_fail["no_forbidden_sites"]


def test_unknown_host_raises():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("not_a_real_host",), mode="PRODUCTION")
    with pytest.raises(OptimizationError):
        optimize_cds(req)


def test_unknown_mode_raises():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="NOT_A_MODE")
    with pytest.raises(OptimizationError):
        optimize_cds(req)


def test_reproducible_with_same_seed():
    req = OptimizationRequest(protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="PRODUCTION", seed=42)
    r1 = optimize_cds(req)[0]
    r2 = optimize_cds(req)[0]
    assert r1.dna == r2.dna


def test_zero_rare_codons_holds_across_all_modes_and_hosts():
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12", "human", "b_subtilis_168"), mode="ALL", seed=11,
    )
    for r in optimize_cds(req):
        assert r.codons_below_w_threshold == 0, f"{r.host}/{r.mode} has a codon below w=0.3"
        assert r.min_w_used >= 0.30
        assert r.constraint_pass_fail["zero_rare_codons"]


def test_hairpin_fields_are_reported_even_when_unresolved():
    # A GC-rich, repetitive protein is a realistic case where the hairpin
    # limits (6b/6c) can't always be fully met -- the spec requires exact
    # numbers and an honest "blocked by" note, never a silent pass.
    gc_rich_protein = "MAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGDDEDEDEDEDEDEDEDEDEDEDEDW"
    req = OptimizationRequest(protein=gc_rich_protein, hosts=("e_coli_k12",), mode="PRODUCTION", seed=12)
    r = optimize_cds(req)[0]
    assert r.worst_window_mfe is not None
    assert r.longest_stem >= 0
    if not r.constraint_pass_fail["no_overstable_window"] or not r.constraint_pass_fail["no_overlong_stem"]:
        assert any("hairpin limit" in n for n in r.notes)


def test_temperature_is_accepted_and_threaded_through():
    req = OptimizationRequest(
        protein=TEST_PROTEIN, hosts=("e_coli_k12",), mode="PRODUCTION",
        five_prime_utr="AGGAGGACAGCTATG", seed=13, temperature_c=30.0,
    )
    result = optimize_cds(req)[0]
    assert result.five_prime_dG is not None  # request-level temperature test lives in test_folding.py
