import pytest

from plasmid_design.codon_usage import (
    DEFAULT_ORGANISM,
    CodonUsageError,
    load_codon_table,
)


def test_loads_default_organism():
    table = load_codon_table(DEFAULT_ORGANISM)
    assert table.organism == "e_coli_k12"
    assert len(table.by_codon) == 64


def test_all_20_amino_acids_and_stop_present():
    table = load_codon_table(DEFAULT_ORGANISM)
    expected = set("ACDEFGHIKLMNPQRSTVWY") | {"*"}
    assert set(table.by_amino_acid) == expected


def test_fractions_sum_to_one_per_amino_acid():
    table = load_codon_table(DEFAULT_ORGANISM)
    for aa, options in table.by_amino_acid.items():
        total = sum(opt.fraction for opt in options)
        assert total == pytest.approx(1.0, abs=0.02), aa


def test_methionine_and_tryptophan_are_unambiguous():
    table = load_codon_table(DEFAULT_ORGANISM)
    (met,) = table.codons_for("M")
    assert met.codon == "ATG"
    assert met.fraction == 1.0
    (trp,) = table.codons_for("W")
    assert trp.codon == "TGG"
    assert trp.fraction == 1.0


def test_unknown_organism_raises():
    with pytest.raises(CodonUsageError):
        load_codon_table("dragon")
