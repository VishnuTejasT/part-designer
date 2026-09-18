import pytest

from plasmid_design.codon_usage import load_codon_table
from plasmid_design.reverse_translate import (
    ProteinSequenceError,
    clean_protein_sequence,
    reverse_translate,
)

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"


def test_clean_protein_sequence_strips_whitespace_and_uppercases():
    assert clean_protein_sequence(" mkt\nAY \t ") == "MKTAY"


def test_reverse_translate_appends_stop_codon_when_absent():
    result = reverse_translate("MKT", seed=1)
    assert result.stop_codon_appended is True
    assert len(result.dna) == 4 * 3  # M K T + stop
    assert result.choices[-1].amino_acid == "*"


def test_reverse_translate_respects_explicit_stop():
    result = reverse_translate("MKT*", seed=1)
    assert result.stop_codon_appended is False
    assert len(result.dna) == 4 * 3


def test_reverse_translate_deterministic_with_seed():
    a = reverse_translate(PROTEIN, seed=42)
    b = reverse_translate(PROTEIN, seed=42)
    assert a.dna == b.dna


def test_reverse_translate_uses_real_codon_frequencies_not_always_the_top_codon():
    # Leucine's most common E. coli codon (CTG, fraction 0.46) should not be
    # chosen every single time across many draws -- that would indicate the
    # implementation is just hardcoding the top codon instead of sampling
    # from the real frequency distribution.
    result = reverse_translate("L" * 200, seed=7)
    leu_codons = [c.codon for c in result.choices if c.amino_acid == "L"]
    assert len(set(leu_codons)) > 1


def test_reverse_translate_round_trips_through_standard_genetic_code():
    table = load_codon_table("e_coli_k12")
    result = reverse_translate(PROTEIN, seed=3)
    translated_back = "".join(
        table.by_codon[result.dna[i : i + 3]]["amino_acid"]
        for i in range(0, len(result.dna), 3)
    )
    assert translated_back == PROTEIN + "*"


def test_rejects_non_amino_acid_characters():
    with pytest.raises(ProteinSequenceError):
        reverse_translate("MKT123")


def test_rejects_internal_stop_codon():
    with pytest.raises(ProteinSequenceError):
        reverse_translate("MK*T")


def test_rejects_empty_sequence():
    with pytest.raises(ProteinSequenceError):
        reverse_translate("")
