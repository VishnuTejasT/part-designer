from plasmid_design.rfc10 import (
    find_site_occurrences,
    reverse_complement,
    validate_rfc10,
)


def test_reverse_complement():
    assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI is palindromic
    assert reverse_complement("AAAACCC") == "GGGTTTT"


def test_clean_sequence_passes_with_no_occurrences():
    # No illegal sites, no NotI site.
    dna = "ATGAAAACCGGATTTAAACGTCTGCCGTAA"
    report = validate_rfc10(dna)
    assert report.passed
    assert report.hard_failures == ()
    assert report.warnings == ()


def test_ecori_forward_site_detected_at_exact_position():
    dna = "ATG" + "GAATTC" + "TAA"
    report = validate_rfc10(dna)
    assert not report.passed
    ecori = next(c for c in report.checks if c.enzyme == "EcoRI")
    assert not ecori.passed
    assert len(ecori.occurrences) == 1
    occ = ecori.occurrences[0]
    assert (occ.start, occ.end) == (4, 9)
    assert occ.matched_sequence == "GAATTC"


def test_all_four_hard_fail_sites_individually_detected():
    hard_fail_sites = {
        "EcoRI": "GAATTC",
        "XbaI": "TCTAGA",
        "SpeI": "ACTAGT",
        "PstI": "CTGCAG",
    }
    for enzyme, site in hard_fail_sites.items():
        dna = "ATGAAA" + site + "AAATAA"
        report = validate_rfc10(dna)
        check = next(c for c in report.checks if c.enzyme == enzyme)
        assert not check.passed, enzyme
        assert not report.passed, enzyme


def test_noti_is_warning_not_hard_fail():
    dna = "ATGAAA" + "GCGGCCGC" + "AAATAA"
    report = validate_rfc10(dna)
    notI = next(c for c in report.checks if c.enzyme == "NotI")
    assert not notI.passed
    assert notI.severity == "warning"
    # NotI alone must not flip the overall pass/fail result.
    assert report.passed
    assert report.hard_failures == ()
    assert len(report.warnings) == 1


def test_reverse_strand_detection_for_non_palindromic_pattern():
    # GATTACA's reverse complement is TGTAATC; embed the reverse complement
    # in the forward sequence and confirm the scanner reports it as a '-'
    # strand hit mapped back to correct forward-strand coordinates.
    pattern = "GATTACA"
    embedded_rc = reverse_complement(pattern)  # TGTAATC
    seq = "AAA" + embedded_rc + "TTT"  # positions 4-10 (1-indexed)
    occurrences = find_site_occurrences(seq, pattern)
    assert len(occurrences) == 1
    occ = occurrences[0]
    assert occ["strand"] == "-"
    assert (occ["start"], occ["end"]) == (4, 10)
    assert occ["matched_sequence"] == embedded_rc


def test_palindromic_hit_reported_once_not_duplicated():
    dna = "ATG" + "GAATTC" + "TAA"
    occurrences = find_site_occurrences(dna, "GAATTC")
    assert len(occurrences) == 1
    assert occurrences[0]["strand"] == "+/-"


def test_multiple_matches_all_reported():
    dna = "CTGCAG" + "CTGCAG"
    occurrences = find_site_occurrences(dna, "CTGCAG")
    assert [(o["start"], o["end"]) for o in occurrences] == [(1, 6), (7, 12)]
