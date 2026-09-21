"""Part Finder. Data checks run against the real extract built from the checksummed
Registry dump; rules are pinned with small synthetic parts that are labeled as such."""
import pytest

from plasmid_design import part_finder as pf
from plasmid_design.part_finder import PartFinderError, find_parts, level_hint, registry_url, rfc10_status

DUMP_SHA256 = "c64bbc9a19b0f84eb314026a5ba509d76239634c5d5061b158f4e555e2fdcd50"
ILLEGAL = ("GAATTC", "TCTAGA", "ACTAGT", "CTGCAG")


@pytest.fixture(scope="module")
def data():
    return pf.load_dataset()


def by_name(data):
    return {p["name"]: p for p in data["parts"]}


# ---- provenance ---------------------------------------------------------------------------
def test_dataset_records_the_exact_dump_it_came_from(data):
    assert data["source"]["sha256"] == DUMP_SHA256
    assert data["source"]["file"].startswith("sql_parts")
    assert data["count"] == len(data["parts"]) > 3000
    assert {p["kind"] for p in data["parts"]} == {"promoter", "rbs", "terminator"}


def test_b0034_matches_the_registry(data):
    p = by_name(data)["BBa_B0034"]
    assert p["sequence"] == "AAAGAGGAGAAA" and p["kind"] == "rbs" and p["status"] == "Available"


def test_every_returned_part_traces_to_the_dataset_and_is_rfc10_clean(data):
    names = by_name(data)
    r = find_parts("e_coli_k12", top=10)
    total = 0
    for kind, block in r["results"].items():
        for part in block["parts"]:
            total += 1
            src = names[part["name"]]                      # traceable
            assert src["kind"] == kind
            assert not any(site in src["sequence"] or site in rc(src["sequence"]) for site in ILLEGAL)
            assert part["registry_url"].startswith("https://registry.igem.org/parts/bba-")
    assert total >= 25


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---- RFC10 filter reuses Part 1's scanner -----------------------------------------------------
def test_parts_with_illegal_sites_are_excluded(data):
    bad = [p for p in data["parts"] if p["kind"] == "promoter" and any(s in p["sequence"] for s in ILLEGAL)]
    assert bad, "dataset should contain some RFC10-illegal promoters for this test to mean anything"
    listed = {x["name"] for x in find_parts("e_coli_k12", kinds=("promoter",), top=10)["results"]["promoter"]["parts"]}
    assert not listed & {p["name"] for p in bad}
    block = find_parts("e_coli_k12", kinds=("promoter",))["results"]["promoter"]
    assert block["passing_filters"] < block["considered"]


def test_rfc10_status_matches_part1_scanner_and_flags_noti():
    assert rfc10_status("ATGGAATTCTAA")["illegal_sites"] == ["EcoRI"]
    assert rfc10_status("ATGCTGCAGTAA")["ok"] is False
    noti = rfc10_status("AAAGCGGCCGCAAA")
    assert noti["ok"] is True and noti["warnings"] == ["NotI"]        # flagged, not blocked
    assert rfc10_status("AAAGAGGAGAAA")["ok"] is True


def test_type_iis_avoidance_is_opt_in():
    seq = "AAAGGTCTCAAA"  # BsaI site: legal for RFC10
    assert rfc10_status(seq)["ok"] is True
    assert rfc10_status(seq, avoid_type_iis=True)["ok"] is False


def test_check_cds_reuses_rfc10():
    assert pf.check_cds("ATGGAATTCTAA") == {"ok": False, "illegal_sites": ["EcoRI"], "warnings": []}
    assert pf.check_cds("ATGAAATAA")["ok"] is True


# ---- chassis ------------------------------------------------------------------------------------
def test_ecoli_results_are_tagged_for_ecoli_or_untagged_never_other_chassis(data):
    names = by_name(data)
    for block in find_parts("e_coli_bl21_de3", top=10)["results"].values():
        for part in block["parts"]:
            tags = names[part["name"]]["chassis"]
            if part["chassis_evidence"] == "tagged":
                assert any(t.startswith("//chassis/prokaryote/ecoli") for t in tags)
            else:
                assert part["chassis_evidence"] == "not_recorded" and tags == []


def test_yeast_only_parts_are_not_offered_for_ecoli(data):
    yeast_only = {p["name"] for p in data["parts"] if p["chassis"] and all("yeast" in t for t in p["chassis"])}
    assert yeast_only
    got = {x["name"] for k in find_parts("e_coli_k12", top=10)["results"].values() for x in k["parts"]}
    assert not got & yeast_only


def test_host_without_a_mapped_tag_is_refused_not_guessed():
    with pytest.raises(PartFinderError, match="No Registry chassis tag"):
        find_parts("c_reinhardtii")


# ---- no invented strength ---------------------------------------------------------------------------
def test_no_part_gets_a_numeric_estimate_without_a_verified_value():
    for block in find_parts("e_coli_k12", top=10)["results"].values():
        for part in block["parts"]:
            e = part["expression"]
            assert e["estimate"] is None
            if e["verified_strength"] is None:
                assert "No numeric expression estimate" in e["note"]


def test_the_one_verified_value_is_exposed_with_its_unit_and_source():
    v = pf.load_verified_strengths()
    assert set(v) == {"BBa_J23100"}
    assert v["BBa_J23100"].efficiency_value == pytest.approx(8.5326)
    part = by_name(pf.load_dataset())["BBa_J23100"]
    c = pf._candidate(part, pf.HOST_CHASSIS["e_coli_k12"], None, v, False)
    assert c["verified"] is not None
    assert "Verified strength value: 8.5326 molecule s^(-1) (per cell)" in " ".join(pf._explain(c, None))


# ---- level hints and warnings ------------------------------------------------------------------------
@pytest.mark.parametrize("text,expected", [
    ("strong constitutive promoter", "high"), ("Weak RBS", "low"), ("medium strength promoter", "moderate"),
    ("strong or weak, depending", None), ("constitutive promoter family member", None), ("", None),
])
def test_level_hint_only_when_the_description_says_so(text, expected):
    h = level_hint(text)
    assert (h["level"] if h else None) == expected
    if h:
        assert h["source"] == "part description text"


def test_strong_plus_strong_warning_is_labeled_as_text_derived():
    hi = {"expression": {"level_hint": {"level": "high"}}}
    lo = {"expression": {"level_hint": {"level": "low"}}}
    w = pf._warnings({"promoter": {"parts": [hi]}, "rbs": {"parts": [hi]}})
    assert len(w) == 1 and "description text" in w[0]
    assert pf._warnings({"promoter": {"parts": [hi]}, "rbs": {"parts": [lo]}}) == []
    assert pf._warnings({"promoter": {"parts": []}, "rbs": {"parts": [hi]}}) == []


def test_level_target_ranks_stated_matches_first_within_availability():
    r = find_parts("e_coli_k12", level="high", kinds=("promoter",), top=10)["results"]["promoter"]["parts"]
    statuses = [p["status"] for p in r]
    assert statuses == sorted(statuses, key=lambda s: s != "Available")      # Available strictly first
    avail = [p for p in r if p["status"] == "Available"]
    hits = [any("matches your target" in x for x in p["reasons"]) for p in avail]
    assert hits == sorted(hits, reverse=True)                                  # matches before non-matches


# ---- output contract -------------------------------------------------------------------------------------
def test_top_n_and_reasons_and_links():
    r = find_parts("e_coli_k12", top=3)
    for block in r["results"].values():
        assert len(block["parts"]) <= 3
        for p in block["parts"]:
            assert p["reasons"] and p["reasons"][0].startswith("Passes the RFC10 check")
            assert p["registry_url"] == registry_url(p["name"])
    assert r["dataset"]["source"]["sha256"] == DUMP_SHA256 and r["limitations"]


def test_registry_url_slug_rule():
    assert registry_url("BBa_J23100") == "https://registry.igem.org/parts/bba-j23100"


@pytest.mark.parametrize("kwargs", [{"standard": "RFC12"}, {"level": "extreme"}, {"top": 0}, {"top": 11}, {"kinds": ("cds",)}])
def test_bad_input_is_rejected(kwargs):
    with pytest.raises(PartFinderError):
        find_parts("e_coli_k12", **kwargs)


def test_regulation_comes_from_registry_tags_and_orders_always_on_first(data):
    names = by_name(data)
    assert pf.regulation_kind(names["BBa_J23100"]) == "constitutive"
    assert pf.regulation_kind(names["BBa_R0040"]) == "regulated"       # //regulation/negative (TetR repressible)
    assert pf.regulation_kind(names["BBa_R0010"]) == "regulated"
    assert pf.regulation_kind(names["BBa_B0034"]) is None               # not a promoter
    assert pf.regulation_kind({"kind": "promoter", "regulation": []}) == "unknown"
    assert pf.regulation_kind({"kind": "promoter", "regulation": ["//regulation/constitutive", "//regulation/negative"]}) == "regulated"
    rank = {"constitutive": 0, "unknown": 1, "regulated": 2}
    groups = {}
    for p in find_parts("e_coli_k12", kinds=("promoter",), top=10)["results"]["promoter"]["parts"]:
        key = (p["status"] == "Available", any("matches your target" in r for r in p["reasons"]))
        groups.setdefault(key, []).append(rank[p["regulation"]])
    assert all(v == sorted(v) for v in groups.values())                 # within a group, always-on comes first


def test_only_promoters_carry_a_regulation_label():
    r = find_parts("e_coli_k12", top=3)["results"]
    assert all(p["regulation"] in ("constitutive", "regulated", "unknown") for p in r["promoter"]["parts"])
    assert all(p["regulation"] is None for k in ("rbs", "terminator") for p in r[k]["parts"])


def test_kinds_subset():
    assert set(find_parts("e_coli_k12", kinds=("rbs",))["results"]) == {"rbs"}


def test_rbs_options_are_real_short_available_ecoli_parts(data):
    opts = pf.rbs_options("e_coli_bl21_de3")
    names = by_name(data)
    assert opts and len(opts) <= 8
    for o in opts:
        src = names[o["name"]]
        assert o["dna"] == src["sequence"] and len(o["dna"]) <= 30 and src["status"] == "Available"
        assert any("ecoli" in t for t in src["chassis"])
    b0034 = next(o for o in opts if o["name"] == "BBa_B0034")
    assert b0034["dna"] == "AAAGAGGAGAAA"                       # the value that used to be hardcoded


def test_find_parts_exposes_sequence_for_verification(data):
    for block in find_parts("e_coli_k12", top=3)["results"].values():
        for p in block["parts"]:
            assert p["sequence"] == by_name(data)[p["name"]]["sequence"]


def test_t7_promoters_carry_the_polymerase_caution():
    part = next(p for p in pf.load_dataset()["parts"] if p["kind"] == "promoter" and "T7" in p["short_desc"]
                and not any(s in p["sequence"] for s in ILLEGAL))
    c = pf._candidate(part, None, None, {}, False)
    assert any("T7 RNA polymerase" in r for r in pf._explain(c, None))
    other = next(p for p in pf.load_dataset()["parts"] if p["kind"] == "promoter" and "T7" not in p["short_desc"]
                 and not any(s in p["sequence"] for s in ILLEGAL))
    assert not any("T7" in r for r in pf._explain(pf._candidate(other, None, None, {}, False), None))
