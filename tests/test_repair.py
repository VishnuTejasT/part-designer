import dataclasses
import random

from plasmid_design import repair as repair_module
from plasmid_design.codon_usage import CodonOption, load_codon_table
from plasmid_design.design import design_cds
from plasmid_design.repair import repair_illegal_sites
from plasmid_design.reverse_translate import CodonChoice, reverse_translate
from plasmid_design.rfc10 import validate_rfc10

# GAA (Glu) + TTC (Phe) spells the EcoRI site GAATTC outright, and Glu/Phe
# both have multiple synonymous codons, so this sequence reliably produces
# a hard-fail site early (seed 0) that also has room to be repaired.
REPEATED_PROTEIN = "EF" * 30


def _first_seed_with_a_hard_fail(protein: str, limit: int = 500) -> int:
    for seed in range(limit):
        if validate_rfc10(reverse_translate(protein, seed=seed).dna).hard_failures:
            return seed
    raise AssertionError(f"no hard-fail site found in {limit} seeds")


def test_repair_removes_illegal_site_without_changing_protein():
    seed = _first_seed_with_a_hard_fail(REPEATED_PROTEIN)

    unrepaired = design_cds(REPEATED_PROTEIN, seed=seed, auto_fix=False)
    repaired = design_cds(REPEATED_PROTEIN, seed=seed, auto_fix=True)

    assert not unrepaired.validation.passed
    assert unrepaired.repair is None

    assert repaired.validation.passed
    assert repaired.repair is not None
    assert repaired.repair.resolved is True
    assert len(repaired.repair.changes) > 0

    # Same protein in, same protein out -- only synonymous codons changed.
    assert repaired.translation.protein == unrepaired.translation.protein
    table = load_codon_table(repaired.translation.organism)
    translated_back = "".join(
        table.by_codon[repaired.translation.dna[i : i + 3]]["amino_acid"]
        for i in range(0, len(repaired.translation.dna), 3)
    )
    assert translated_back == repaired.translation.protein + "*"


def test_repair_result_resolved_flag_matches_final_validation():
    seed = _first_seed_with_a_hard_fail(REPEATED_PROTEIN)
    result = design_cds(REPEATED_PROTEIN, seed=seed, auto_fix=True)
    assert result.repair.resolved == result.validation.passed


def test_repair_is_noop_when_sequence_already_clean():
    result = design_cds("MKT", seed=1, auto_fix=True)
    assert result.validation.passed
    assert result.repair is None


def test_repair_gives_up_gracefully_when_no_synonymous_alternative_exists(monkeypatch):
    """Stub the codon table so the only 'option' for each amino acid is the
    codon already in place -- i.e. no synonymous change is possible -- and
    confirm repair leaves the sequence untouched and reports unresolved
    instead of looping or raising."""

    translation = reverse_translate("MW", seed=1)
    forced_choices = (
        CodonChoice(position=1, amino_acid="M", codon="GAA", fraction=1.0),
        CodonChoice(position=2, amino_acid="W", codon="TTC", fraction=1.0),
    )
    broken = dataclasses.replace(translation, dna="GAATTC", choices=forced_choices)

    class StubTable:
        _codon = {"M": "GAA", "W": "TTC"}

        def codons_for(self, amino_acid):
            return (CodonOption(codon=self._codon[amino_acid], fraction=1.0, per_thousand=0.0),)

    monkeypatch.setattr(repair_module, "load_codon_table", lambda organism: StubTable())

    result = repair_illegal_sites(broken, random.Random(1))

    assert result.resolved is False
    assert result.changes == ()
    assert result.translation.dna == "GAATTC"
    assert result.iterations == 1


def test_codon_positions_overlapping():
    overlap = repair_module._codon_positions_overlapping
    assert overlap(1, 6) == [1, 2]  # exactly two codons, no offset
    assert overlap(4, 9) == [2, 3]
    assert overlap(3, 8) == [1, 2, 3]  # offset by one base, spans three codons
