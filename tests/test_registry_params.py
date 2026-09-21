"""Efficiency parser. The J23100 test uses a fixture transcribed from the real
Registry XML (see the fixture's header comment). Every other case is a
SYNTHETIC input built here to pin one rule; none of it is part data."""
from pathlib import Path

import pytest

from plasmid_design.registry_params import (
    ABSENT, NA_ONLY, VALID, coverage, load_xml_dir, normalize_unit, parse_efficiency_value, parse_part_xml,
)

FIX = Path(__file__).parent / "fixtures" / "registry"


def part_xml(name="BBa_SYN", rows=(), extra=""):
    """One synthetic <part>. rows: (name, value, id, m_date[, units])."""
    ps = "".join(
        f"<parameter><name>{r[0]}</name><value>{r[1]}</value><units>{r[4] if len(r) > 4 else ''}</units>"
        f"<id>{r[2]}</id><m_date>{r[3]}</m_date></parameter>" for r in rows
    )
    return (f"<part><part_id>1</part_id><part_name>{name}</part_name><part_type>Regulatory</part_type>"
            f"<sequences><seq_data>acgt</seq_data></sequences><parameters>{ps}</parameters>{extra}</part>")


def wrap(*parts):
    return "<rsbpml><part_list>" + "".join(parts) + "</part_list></rsbpml>"


def one(*rows, **kw):
    return parse_part_xml(wrap(part_xml(rows=rows, **kw)))[0]


# ---- the real J23100 structure ------------------------------------------------------
def test_real_j23100_fixture_is_valid_with_the_newer_numeric_row():
    (rec,) = parse_part_xml((FIX / "BBa_J23100_transcribed_from_pdf.xml").read_text())
    assert rec.part_name == "BBa_J23100" and rec.part_id == 7142
    assert rec.efficiency_status == VALID
    assert rec.efficiency_value == pytest.approx(8.5326)
    assert rec.efficiency_unit == "molecule s^(-1) (per cell)"
    assert rec.efficiency.row_id == 30083          # newer than the NA row 30082
    assert len(rec.efficiency_rows) == 2           # both efficiency rows kept; other names ignored
    assert rec.newer_non_numeric_after_valid is False
    assert rec.sequence == "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC"
    assert (rec.release_status, rec.sample_status, rec.part_results) == ("Released HQ 2013", "In stock", "Works")


# ---- value parsing ---------------------------------------------------------------------
@pytest.mark.parametrize("raw,number,unit", [
    ("8.5326 molecule s^(−1 ) (per cell)", 8.5326, "molecule s^(-1) (per cell)"),
    ("0.5", 0.5, ""),
    ("  12  au ", 12.0, "au"),
    ("1.2e3 REU", 1200.0, "REU"),
    ("-0.4 x", -0.4, "x"),
    ("−0.4 x", -0.4, "x"),
])
def test_numeric_values_parse(raw, number, unit):
    assert parse_efficiency_value(raw) == (number, unit)


@pytest.mark.parametrize("raw", ["NA", "n/a", "", "   ", "high", "strong promoter 5", "approx 3", "e5", None])
def test_non_numeric_values_are_rejected(raw):
    assert parse_efficiency_value(raw) == (None, "")


def test_units_element_used_only_when_value_has_no_inline_unit():
    assert parse_efficiency_value("3", "PoPS") == (3.0, "PoPS")
    assert parse_efficiency_value("3 REU", "PoPS") == (3.0, "REU")


def test_unit_spellings_normalize_to_the_same_string():
    assert normalize_unit("molecule s^(−1 ) (per cell)") == normalize_unit("molecule  s^(-1)  (per cell)")


# ---- filtering the noisy history ---------------------------------------------------------
def test_only_the_exact_name_efficiency_counts():
    rec = one(("Efficiency", "5", 1, "2020-01-01 00:00:00"), ("efficiency_note", "5", 2, "2020-01-01 00:00:00"),
              ("n/a", "5", 3, "2020-01-01 00:00:00"))
    assert rec.efficiency_status == ABSENT and rec.efficiency_rows == ()


def test_no_efficiency_row_is_absent():
    assert one(("positive_regulators", "BBa_B0034", 1, "2009-01-01 00:00:00")).efficiency_status == ABSENT


def test_only_na_rows_is_na_only():
    rec = one(("efficiency", "NA", 1, "2020-01-01 00:00:00"), ("efficiency", "n/a", 2, "2021-01-01 00:00:00"))
    assert rec.efficiency_status == NA_ONLY and rec.efficiency is None and rec.efficiency_value is None


def test_newest_numeric_row_wins_by_date_then_id():
    rec = one(("efficiency", "1 x", 5, "2019-01-01 00:00:00"), ("efficiency", "2 x", 9, "2021-01-01 00:00:00"),
              ("efficiency", "3 x", 7, "2020-01-01 00:00:00"))
    assert rec.efficiency_value == 2.0
    tie = one(("efficiency", "1 x", 5, "2020-01-01 00:00:00"), ("efficiency", "2 x", 6, "2020-01-01 00:00:00"))
    assert tie.efficiency_value == 2.0


def test_a_newer_na_after_a_valid_value_is_flagged_not_hidden():
    rec = one(("efficiency", "4 x", 1, "2019-01-01 00:00:00"), ("efficiency", "NA", 2, "2022-01-01 00:00:00"))
    assert rec.efficiency_status == VALID and rec.efficiency_value == 4.0
    assert rec.newer_non_numeric_after_valid is True


def test_multiple_parts_in_one_part_list_and_dir_loading(tmp_path):
    xml = wrap(part_xml("BBa_A", rows=[("efficiency", "1 x", 1, "2020-01-01 00:00:00")]), part_xml("BBa_B"))
    (tmp_path / "part.BBa_A.xml").write_text(xml)
    (tmp_path / "ignored.txt").write_text("x")
    recs = load_xml_dir(tmp_path)
    assert [r.part_name for r in recs] == ["BBa_A", "BBa_B"]
    assert [r.efficiency_status for r in recs] == [VALID, ABSENT]


def test_malformed_xml_raises_instead_of_silently_returning_nothing():
    import xml.etree.ElementTree as ET
    with pytest.raises(ET.ParseError):
        parse_part_xml("<rsbpml><part_list><part>")


# ---- coverage ---------------------------------------------------------------------------------
def test_coverage_counts_status_kind_and_units():
    recs = parse_part_xml(wrap(
        part_xml("BBa_P1", rows=[("efficiency", "1 REU", 1, "2020-01-01 00:00:00")]),
        part_xml("BBa_P2", rows=[("efficiency", "2 REU", 1, "2020-01-01 00:00:00")]),
        part_xml("BBa_P3", rows=[("efficiency", "3 au", 1, "2020-01-01 00:00:00")]),
        part_xml("BBa_P4", rows=[("efficiency", "NA", 1, "2020-01-01 00:00:00")]),
        part_xml("BBa_R1"),
    ))
    kinds = {"BBa_P1": "promoter", "BBa_P2": "promoter", "BBa_P3": "promoter", "BBa_P4": "promoter", "BBa_R1": "rbs"}
    c = coverage(recs, kinds)
    assert c["total"] == 5 and c["status"] == {VALID: 3, NA_ONLY: 1, ABSENT: 1}
    assert c["by_kind"]["promoter"] == {VALID: 3, NA_ONLY: 1, ABSENT: 0}
    assert c["by_kind"]["rbs"] == {VALID: 0, NA_ONLY: 0, ABSENT: 1}
    assert c["kind_units"]["promoter"] == {"REU": 2, "au": 1}   # rankable only within one unit
