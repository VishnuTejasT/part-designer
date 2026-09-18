import pytest

from plasmid_design.codon_usage import load_codon_table
from plasmid_design.trna_usage import (
    TRNA_TABLES,
    TrnaDataError,
    load_trna_table,
    relative_adaptiveness,
    tai_score,
)


def test_trna_tables_registry_matches_confirmed_hosts():
    assert set(TRNA_TABLES) == {"e_coli_k12", "s_cerevisiae", "b_subtilis_168", "human"}


def test_unavailable_host_raises_clear_error():
    with pytest.raises(TrnaDataError):
        load_trna_table("synechocystis_pcc6803")


def test_e_coli_gene_count_matches_known_value():
    # 87 tRNA genes in E. coli K-12 is a well-known, independently-cited
    # number (also the exact total in dos Reis' own reference dataset).
    table = load_trna_table("e_coli_k12")
    assert sum(table.anticodon_counts.values()) == 87


def test_relative_adaptiveness_scores_between_zero_and_one():
    trna_table = load_trna_table("e_coli_k12")
    codon_table = load_codon_table("e_coli_k12")
    w = relative_adaptiveness(trna_table, codon_table)
    assert len(w) == 60  # 64 codons minus 3 stops minus Met
    assert all(0 < v <= 1.0 for v in w.values())
    assert max(w.values()) == pytest.approx(1.0)


def test_tai_score_of_max_w_codons_is_one():
    trna_table = load_trna_table("e_coli_k12")
    codon_table = load_codon_table("e_coli_k12")
    w = relative_adaptiveness(trna_table, codon_table)
    best_codon = max(w, key=w.get)
    dna = best_codon * 5
    assert tai_score(dna, trna_table, codon_table) == pytest.approx(1.0, abs=1e-6)
