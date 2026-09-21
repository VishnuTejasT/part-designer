"""Parse per-part Registry XML (rsbpml) and extract the ``efficiency`` parameter.

The Registry keeps years of edit history in <parameters>: many NA / n/a /
duplicate rows per part. Only rows named exactly ``efficiency`` whose value
starts with a parseable number count. Where several such rows exist the newest
one (by m_date, then id) wins, and a newer non-numeric row is recorded rather
than hidden.

Reads local files only (the Registry blocks scripted downloads); nothing here
touches the network.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

NUMBER_FIRST = re.compile(r"^\s*([-+−]?\d+(?:\.\d+)?(?:[eE][-+−]?\d+)?)\s*(.*?)\s*$", re.S)

VALID, NA_ONLY, ABSENT = "valid", "na_only", "absent"


@dataclass(frozen=True)
class EfficiencyRow:
    row_id: int
    m_date: str
    raw_value: str
    number: float | None
    unit: str  # normalized; "" when the row gives a bare number


@dataclass(frozen=True)
class PartRecord:
    part_name: str
    part_id: int | None
    part_type: str
    release_status: str
    sample_status: str
    part_results: str
    part_url: str
    sequence: str
    efficiency_status: str  # valid | na_only | absent
    efficiency: EfficiencyRow | None  # the newest numeric row, if any
    efficiency_rows: tuple[EfficiencyRow, ...] = field(default=())
    newer_non_numeric_after_valid: bool = False

    @property
    def efficiency_value(self) -> float | None:
        return self.efficiency.number if self.efficiency else None

    @property
    def efficiency_unit(self) -> str:
        return self.efficiency.unit if self.efficiency else ""


def normalize_unit(text: str) -> str:
    """'molecule s^(−1 ) (per cell)' -> 'molecule s^(-1) (per cell)'. Same
    unit written slightly differently must compare equal."""
    t = text.replace("−", "-").replace("–", "-")
    t = re.sub(r"\s+", " ", t).strip()
    t = re.sub(r"\(\s+", "(", t)
    t = re.sub(r"\s+\)", ")", t)
    return t


def parse_efficiency_value(raw: str, units_element: str = "") -> tuple[float | None, str]:
    """(number, unit) for a value such as '8.5326 molecule s^(-1) (per cell)';
    (None, '') for NA, n/a, empty, or anything that doesn't start with a number."""
    if raw is None:
        return None, ""
    m = NUMBER_FIRST.match(raw)
    if not m:
        return None, ""
    try:
        number = float(m.group(1).replace("−", "-"))
    except ValueError:
        return None, ""
    unit = normalize_unit(m.group(2)) or normalize_unit(units_element or "")
    return number, unit


def _text(node: ET.Element | None, tag: str) -> str:
    child = node.find(tag) if node is not None else None
    return (child.text or "").strip() if child is not None else ""


def _rows(part: ET.Element) -> list[EfficiencyRow]:
    rows = []
    for p in part.findall("./parameters/parameter"):
        if _text(p, "name") != "efficiency":
            continue
        raw = _text(p, "value")
        number, unit = parse_efficiency_value(raw, _text(p, "units"))
        try:
            rid = int(_text(p, "id"))
        except ValueError:
            rid = -1
        rows.append(EfficiencyRow(rid, _text(p, "m_date"), raw, number, unit))
    return rows


def _record(part: ET.Element) -> PartRecord:
    rows = _rows(part)
    ordered = sorted(rows, key=lambda r: (r.m_date, r.row_id), reverse=True)  # newest first
    numeric = [r for r in ordered if r.number is not None]
    best = numeric[0] if numeric else None
    newer_bad = bool(best and any((r.m_date, r.row_id) > (best.m_date, best.row_id) and r.number is None for r in rows))
    status = VALID if best else (NA_ONLY if rows else ABSENT)
    pid = _text(part, "part_id")
    seq = "".join(_text(part, "sequences/seq_data").split()).upper()
    return PartRecord(
        part_name=_text(part, "part_name"), part_id=int(pid) if pid.isdigit() else None,
        part_type=_text(part, "part_type"), release_status=_text(part, "release_status"),
        sample_status=_text(part, "sample_status"), part_results=_text(part, "part_results"),
        part_url=_text(part, "part_url"), sequence=seq, efficiency_status=status, efficiency=best,
        efficiency_rows=tuple(rows), newer_non_numeric_after_valid=newer_bad,
    )


def parse_part_xml(xml_text: str) -> list[PartRecord]:
    """All <part> elements in an rsbpml document (one part per file, or a
    <part_list> holding many)."""
    root = ET.fromstring(xml_text)
    return [_record(p) for p in root.iter("part")]


def load_xml_dir(directory: str | Path) -> list[PartRecord]:
    records: list[PartRecord] = []
    for path in sorted(Path(directory).glob("*.xml")):
        records.extend(parse_part_xml(path.read_text(encoding="utf-8")))
    return records


def coverage(records: list[PartRecord], kind_of: dict[str, str] | None = None) -> dict:
    """Counts of valid / na_only / absent overall, per kind (if ``kind_of``
    maps part_name -> promoter|rbs|terminator), and the units seen. Ranking
    is only meaningful within one (kind, unit) group, so those counts are
    reported too."""
    out: dict = {"total": len(records), "status": {}, "by_kind": {}, "units": {}, "kind_units": {}}
    for r in records:
        out["status"][r.efficiency_status] = out["status"].get(r.efficiency_status, 0) + 1
        kind = (kind_of or {}).get(r.part_name, "unknown")
        bk = out["by_kind"].setdefault(kind, {VALID: 0, NA_ONLY: 0, ABSENT: 0})
        bk[r.efficiency_status] += 1
        if r.efficiency_status == VALID:
            out["units"][r.efficiency_unit] = out["units"].get(r.efficiency_unit, 0) + 1
            ku = out["kind_units"].setdefault(kind, {})
            ku[r.efficiency_unit] = ku.get(r.efficiency_unit, 0) + 1
    return out
