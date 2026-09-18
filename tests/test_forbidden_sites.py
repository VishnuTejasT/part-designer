from plasmid_design.forbidden_sites import scan_forbidden_sites


def test_always_on_type_iis_detected():
    dna = "ATG" + "GGTCTC" + "TAA"  # BsaI
    hits = scan_forbidden_sites(dna)
    assert "BsaI" in hits


def test_rfc10_hard_fail_site_also_detected():
    dna = "ATG" + "GAATTC" + "TAA"  # EcoRI
    hits = scan_forbidden_sites(dna)
    assert "EcoRI" in hits


def test_rfc10_warning_only_site_not_included():
    dna = "ATG" + "GCGGCCGC" + "TAA"  # NotI, warning-severity in RFC10
    hits = scan_forbidden_sites(dna)
    assert "NotI" not in hits


def test_extra_caller_supplied_site():
    dna = "ATG" + "GATATC" + "TAA"  # EcoRV, not in any built-in list
    assert scan_forbidden_sites(dna) == {}
    hits = scan_forbidden_sites(dna, extra_sites=(("EcoRV", "GATATC"),))
    assert "EcoRV" in hits


def test_clean_sequence_has_no_hits():
    dna = "ATGAAACCCGGGTTTTAA"
    assert scan_forbidden_sites(dna) == {}
