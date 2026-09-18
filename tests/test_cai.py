import pytest

from plasmid_design.cai import cai_score, expression_weighted_cai_score, w_scores
from plasmid_design.codon_usage import load_codon_table


def test_w_scores_max_codon_is_one():
    table = load_codon_table("e_coli_k12")
    scores = w_scores(table)
    for amino_acid, options in table.by_amino_acid.items():
        best = max(options, key=lambda o: o.fraction)
        assert scores[best.codon] == pytest.approx(1.0)


def test_cai_score_of_all_max_codons_is_one():
    table = load_codon_table("e_coli_k12")
    # CTG is E. coli's most-used Leu codon, GGC its most-used Gly codon.
    dna = "CTG" * 3 + "GGC" * 3 + "TAA"
    assert cai_score(dna, table) == pytest.approx(1.0, abs=1e-6)


def test_cai_score_excludes_met_and_trp():
    table = load_codon_table("e_coli_k12")
    with pytest.raises(ValueError):
        cai_score("ATG" + "TGG" + "TAA", table)  # Met, Trp, stop -- nothing scorable


def test_expression_weighted_cai_always_unavailable_for_now():
    table = load_codon_table("e_coli_k12")
    value, reason = expression_weighted_cai_score("ATGAAATAA", table)
    assert value is None
    assert "expression-weighted" in reason
