"""Assembler tests. Junction cases use SYNTHETIC parts (patched in) that are each RFC10-clean
on their own, so a hit can only come from the junction. Real-data cases use the Registry extract."""
import pytest

from plasmid_design import vectorizer as vz
from plasmid_design.rfc10 import validate_rfc10
from plasmid_design.vectorizer import VectorizerError, assemble

CDS = "ATGAAAGTTCTGGCGTAA"


def patch_parts(monkeypatch, **seqs):
    real = vz._part
    def fake(name, kind):
        return {"sequence": seqs[name]} if name in seqs else real(name, kind)
    monkeypatch.setattr(vz, "_part", fake)


def real(**kw):
    return assemble("BBa_J23100", "BBa_B0034", kw.pop("cds", CDS), "BBa_B0015", **kw)


def test_scars_are_the_registry_standard():
    assert vz.SCAR_STANDARD == "TACTAGAG" and vz.SCAR_RBS_CDS == "TACTAG"


def test_backbone_is_the_registry_psb1c3_with_provenance():
    b = vz.load_backbones()["pSB1C3"]
    assert b["length"] == len(b["sequence"]) == 2070
    assert b["sequence"].endswith("GAATTCGCGGCCGCTTCTAGAG") and b["sequence"].startswith("TACTAGTAGCGGCCGCTGCAG")
    assert b["source"]["sha256"].startswith("c64bbc9a19b0")


def test_assembly_order_scars_and_coordinates():
    r = real()
    seq, B = r["sequence"], vz.load_backbones()["pSB1C3"]["sequence"]
    assert seq.startswith(B) and seq[len(B):] == r["insert"]
    kinds = [x["kind"] for x in r["layout"]]
    assert kinds == ["promoter", "scar", "rbs", "scar", "cds", "scar", "terminator"]
    scars = [seq[x["start"] - 1:x["end"]] for x in r["layout"] if x["kind"] == "scar"]
    assert scars == ["TACTAGAG", "TACTAG", "TACTAGAG"]                      # promoter|RBS, RBS|CDS, CDS|terminator
    cds = next(x for x in r["layout"] if x["kind"] == "cds")
    assert seq[cds["start"] - 1:cds["end"]] == CDS
    assert r["size"]["total_bp"] == len(seq) == 2070 + len(r["insert"])
    assert r["fasta"].startswith(">construct_pSB1C3") and r["fasta"].replace("\n", "")[len(">construct_pSB1C3_insert"):] == seq


def test_clean_construct_passes_even_though_the_backbone_carries_flank_sites():
    r = real()
    assert validate_rfc10(r["sequence"]).hard_failures                       # a naive whole-plasmid scan always fails
    assert r["junction_check"]["ok"] and r["junction_check"]["violations"] == []


def test_a_site_created_only_at_a_junction_is_caught(monkeypatch):
    patch_parts(monkeypatch, P_SYN="AAAACTAG")                               # ...ACTAG + scar's T = ACTAGT (SpeI)
    assert validate_rfc10("AAAACTAG").passed                                  # clean on its own
    r = assemble("P_SYN", "BBa_B0034", CDS, "BBa_B0015")
    v = r["junction_check"]["violations"]
    assert not r["junction_check"]["ok"] and [(x["enzyme"], x["at_junction"]) for x in v] == [("SpeI", True)]


def test_a_site_across_the_circular_seam_is_caught(monkeypatch):
    patch_parts(monkeypatch, T_SYN="GGGGACTAG")                              # end of insert + backbone start T = ACTAGT
    assert validate_rfc10("GGGGACTAG").passed
    r = assemble("BBa_J23100", "BBa_B0034", CDS, "T_SYN")
    v = r["junction_check"]["violations"]
    assert not r["junction_check"]["ok"] and v[0]["enzyme"] == "SpeI" and v[0]["at_junction"]


def test_an_illegal_site_inside_the_cds_is_caught_and_not_marked_as_junction():
    r = real(cds="ATG" + "GCT" * 10 + "GAATTC" + "GCT" * 10 + "TAA")            # site well away from every scar
    v = r["junction_check"]["violations"]
    assert not r["junction_check"]["ok"] and [x["enzyme"] for x in v] == ["EcoRI"] and v[0]["at_junction"] is False


def test_type_iis_check_is_opt_in():
    cds = "ATGGGTCTCGCGTAA"                                                   # BsaI, legal for RFC10
    assert real(cds=cds)["junction_check"]["ok"]
    r = real(cds=cds, avoid_type_iis=True)
    assert not r["junction_check"]["ok"] and r["junction_check"]["type_iis_sites"] == ["BsaI"]


def test_gc_flags_use_twists_window_rule_and_label_the_project_range_separately(monkeypatch):
    r = real()
    assert (r["gc"]["twist_high_complexity"], r["gc"]["outside_project_range"]) == (False, False)
    patch_parts(monkeypatch, P_GC="G" * 60)
    hi = assemble("P_GC", "BBa_B0034", CDS, "BBa_B0015")["gc"]
    assert hi["twist_high_complexity"] and hi["window_max"] > 90
    assert "not a vendor limit" in hi["project_note"] and "10%" in hi["twist_note"]


def test_size_is_reported_with_evidence_not_an_invented_cutoff():
    s = real()["size"]
    assert s["within_studied_range"] and "trend, not a pass/fail line" in s["evidence"] and "16.1 kb" in s["evidence"]
    big = real(cds="ATG" + "GCT" * 5000 + "TAA")["size"]
    assert not big["within_studied_range"] and "larger than the studied range" in big["evidence"]


def test_insert_beyond_the_clonal_synthesis_limit_is_flagged():
    assert real()["synthesis"]["beyond_clonal_gene_limit"] is False
    assert real(cds="ATG" + "GCT" * 2400 + "TAA")["synthesis"]["beyond_clonal_gene_limit"] is True


def test_ready_to_order_is_a_single_gate_over_the_real_checks():
    r = real()
    assert r["ready_to_order"] == {"ok": True, "blockers": [], "note": "Nothing here blocks ordering this construct as designed."}
    big = real(cds="ATG" + "GCT" * 2400 + "TAA")
    assert big["ready_to_order"]["ok"] is False
    assert any("7,000 bp" in b for b in big["ready_to_order"]["blockers"])


def test_ready_to_order_blocks_on_a_junction_violation(monkeypatch):
    patch_parts(monkeypatch, P_SYN="AAAACTAG")
    r = assemble("P_SYN", "BBa_B0034", CDS, "BBa_B0015")
    assert r["ready_to_order"]["ok"] is False
    assert any("join" in b for b in r["ready_to_order"]["blockers"])


@pytest.mark.parametrize("kwargs,match", [
    ({"cds": "ATGXXTAA"}, "A, C, G, T"), ({"cds": "ATGAAATA"}, "whole number of codons"), ({"cds": "AAAAAATAA"}, "start with ATG"),
    ({"backbone": "pUC19"}, "Unknown backbone"),
])
def test_bad_input_is_rejected(kwargs, match):
    with pytest.raises(VectorizerError, match=match):
        real(**kwargs)


def test_parts_must_be_the_right_kind_and_exist():
    with pytest.raises(VectorizerError, match="is a rbs, not a promoter"):
        assemble("BBa_B0034", "BBa_B0034", CDS, "BBa_B0015")
    with pytest.raises(VectorizerError, match="not in the Registry extract"):
        assemble("BBa_NOPE", "BBa_B0034", CDS, "BBa_B0015")
