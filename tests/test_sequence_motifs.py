from plasmid_design.sequence_motifs import scan_motifs


def test_prokaryote_sd_like_motif_detected():
    dna = "ATG" + "GGAGG" + "AAATAA"
    hits = scan_motifs(dna, "prokaryote")
    assert any(h.kind == "SD-like" for h in hits)


def test_eukaryote_scan_ignores_sd_like_and_catches_eukaryote_motifs():
    dna = "ATG" + "GGAGG" + "AATAAA" + "TAA"
    prok_hits = scan_motifs(dna, "prokaryote")
    euk_hits = scan_motifs(dna, "eukaryote")
    assert any(h.kind == "SD-like" for h in prok_hits)
    assert not any(h.kind == "SD-like" for h in euk_hits)
    assert any(h.kind == "polyA signal" for h in euk_hits)


def test_au_rich_element_detected_in_eukaryote_scan():
    dna = "ATG" + "ATTTA" + "TAA"
    hits = scan_motifs(dna, "eukaryote")
    assert any(h.kind == "AU-rich element" for h in hits)


def test_unknown_domain_returns_no_hits():
    dna = "ATG" + "GGAGG" + "AATAAA" + "TAA"
    assert scan_motifs(dna, None) == []
