import pytest

RNA = pytest.importorskip("RNA", reason="ViennaRNA not installed in this environment")

from plasmid_design import folding


def test_mfe_returns_structure_and_energy():
    result = folding.mfe("GGGAAACCC")
    assert result is not None
    structure, energy = result
    assert len(structure) == 9
    assert energy <= 0


def test_initiation_window_pads_with_utr_and_cds():
    window = folding.initiation_window("ATGAAACTGGTC", "AGGAGGACAG", upstream=5, downstream=6)
    assert window == "GACAG" + "ATGAAA"


def test_initiation_dG_is_a_float_when_available():
    dG = folding.initiation_dG("ATGAAACTGGTCTAA", "AGGAGGACAGCTATG")
    assert isinstance(dG, float)


def test_initiation_dG_none_for_empty_window():
    assert folding.initiation_dG("", "") is None


def test_temperature_changes_mfe_for_a_fixed_sequence():
    # A real hairpin-forming sequence: colder folding must never be less
    # stable (less negative dG) than hotter folding of the *same* sequence.
    seq = "GGGGGGGAAAACCCCCCC"
    _structure_cold, energy_cold = folding.mfe(seq, temperature_c=4.0)
    _structure_hot, energy_hot = folding.mfe(seq, temperature_c=60.0)
    assert energy_cold <= energy_hot


def test_initiation_unpaired_ok_returns_bool_or_none():
    result = folding.initiation_unpaired_ok("ATGAAACTGGTCTAA", "AGGAGGACAGCTATG")
    assert result in (True, False)
    assert folding.initiation_unpaired_ok("", "") is None
