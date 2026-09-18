import dataclasses

from plasmid_design.report import build_report, format_report_text
from plasmid_design.reverse_translate import reverse_translate
from plasmid_design.rfc10 import validate_rfc10


def test_build_report_structure_and_text_rendering():
    translation = reverse_translate("MKT", seed=1)
    validation = validate_rfc10(translation.dna)
    report = build_report(translation, validation)

    assert report["protein"]["sequence"] == "MKT"
    assert report["dna"]["sequence"] == translation.dna
    assert len(report["codon_choices"]) == 4  # M, K, T, stop
    assert report["validation"]["overall_pass"] == validation.passed

    text = format_report_text(report)
    assert "MKT" in text
    assert translation.dna in text
    assert "RFC10 COMPATIBILITY VALIDATION" in text


def test_report_surfaces_illegal_site_position_in_text():
    translation = reverse_translate("MKT", seed=1)
    # Force in an EcoRI site to check the report surfaces its exact position.
    dna_with_site = translation.dna[:3] + "GAATTC" + translation.dna[3:]
    translation = dataclasses.replace(translation, dna=dna_with_site)
    validation = validate_rfc10(dna_with_site)
    report = build_report(translation, validation)
    text = format_report_text(report)
    assert "EcoRI" in text
    assert "position 4-9" in text
