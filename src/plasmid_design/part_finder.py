"""Part Finder: recommend promoter / RBS / terminator parts for a host.

Every part comes from the local Registry extract (data/registry/parts_dataset.json.gz,
built from a checksummed SQL dump) -- nothing is hardcoded. Filtering reuses
Part 1's scanners (rfc10.validate_rfc10, forbidden_sites.scan_forbidden_sites).

There is NO numeric expression prediction. The Registry dump and the current
REST API carry no structured strength data; the only verified value we hold is
the ``efficiency`` of parts whose legacy XML was saved under data/registry/xml/.
Ordering therefore uses evidence we can show (availability, chassis tag,
'Works' result, community use, verified strength if any) and says so.
"""

from __future__ import annotations

import gzip
import json
import math
import re
from functools import lru_cache

from .codon_usage import DATA_DIR
from .forbidden_sites import scan_forbidden_sites
from .registry_params import load_xml_dir
from .rfc10 import validate_rfc10

DATASET_PATH = DATA_DIR / "registry" / "parts_dataset.json.gz"
XML_DIR = DATA_DIR / "registry" / "xml"
KINDS = ("promoter", "rbs", "terminator")
LEVELS = ("low", "moderate", "high")

# Host -> Registry chassis category tag. Only tags that occur in the dump.
HOST_CHASSIS = {
    "e_coli_k12": "//chassis/prokaryote/ecoli",
    "e_coli_bl21_de3": "//chassis/prokaryote/ecoli",
    "s_cerevisiae": "//chassis/eukaryote/yeast",
    "b_subtilis_168": "//chassis/prokaryote/bsubtilis",
    "human": "//chassis/eukaryote/human",
}

_HINTS = (
    ("high", re.compile(r"\b(strong|high)\b", re.I)),
    ("low", re.compile(r"\b(weak|low)\b", re.I)),
    ("moderate", re.compile(r"\b(medium|moderate)\b", re.I)),
)


class PartFinderError(ValueError):
    pass


@lru_cache(maxsize=1)
def load_dataset() -> dict:
    with gzip.open(DATASET_PATH, "rt", encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def load_verified_strengths() -> dict:
    """part name -> PartRecord, for parts with a valid legacy-XML efficiency."""
    if not XML_DIR.exists():
        return {}
    return {r.part_name: r for r in load_xml_dir(XML_DIR) if r.efficiency_status == "valid"}


def registry_url(name: str) -> str:
    """Current Registry page (slug rule checked against the live API: BBa_J23100 -> bba-j23100)."""
    return "https://registry.igem.org/parts/" + name.lower().replace("_", "-")


def level_hint(short_desc: str) -> dict | None:
    """A coarse level ONLY when the description literally says so; never both."""
    found = [(lvl, m.group(0)) for lvl, rx in _HINTS if (m := rx.search(short_desc or ""))]
    if len(found) != 1:
        return None
    return {"level": found[0][0], "word": found[0][1].lower(), "source": "part description text"}


def rfc10_status(sequence: str, avoid_type_iis: bool = False) -> dict:
    report = validate_rfc10(sequence)  # same scanner Part 1 uses
    bad = [c.enzyme for c in report.hard_failures]
    warn = [c.enzyme for c in report.warnings]
    extra = sorted(scan_forbidden_sites(sequence)) if avoid_type_iis else []
    type_iis = [e for e in extra if e in ("BsaI", "BsmBI", "BbsI", "SapI")]
    return {"ok": not bad and not type_iis, "illegal_sites": bad, "warnings": warn, "type_iis_sites": type_iis}


def _candidate(part: dict, chassis_tag: str | None, target: str | None, strengths: dict, avoid_type_iis: bool):
    rfc = rfc10_status(part["sequence"], avoid_type_iis)
    if not rfc["ok"]:
        return None
    tags = part["chassis"]
    if chassis_tag is None:
        evidence = "not_requested"
    elif any(t == chassis_tag or t.startswith(chassis_tag + "/") for t in tags):
        evidence = "tagged"
    elif not tags:
        evidence = "not_recorded"
    else:
        return None  # tagged only for other chassis: no evidence for this one
    hint = level_hint(part["short_desc"])
    verified = strengths.get(part["name"])
    match = bool(target and hint and hint["level"] == target)
    return {"part": part, "rfc": rfc, "evidence": evidence, "hint": hint, "verified": verified, "match": match}


def _sort_key(c: dict):
    """Available first, then a stated-level match, then host tag, Registry
    'results', community use. A verified strength value is shown but does not
    raise rank: one part's number can't be compared to your target level."""
    p = c["part"]
    return (
        -int(p["status"] == "Available"), -int(c["match"]), -int(c["evidence"] == "tagged"),
        {"Works": 0, "": 1, "None": 1, "Issues": 2, "Fails": 3}.get(p["works"], 1),
        -(p["uses"] or 0), p["name"],
    )


def _explain(c: dict, target: str | None) -> list[str]:
    p, out = c["part"], []
    out.append("Passes the RFC10 check: none of EcoRI, XbaI, SpeI, PstI in the sequence"
               + (f" (NotI present: flagged, not blocked)" if "NotI" in c["rfc"]["warnings"] else ""))
    out.append({"tagged": "Registry categories tag it for this host",
                "not_recorded": "No host is recorded for this part, so compatibility is not confirmed",
                "not_requested": "No host filter applied"}[c["evidence"]])
    out.append("Registry status: Available" if p["status"] == "Available" else f"Registry status: {p['status']} (not currently available)")
    if p["works"] in ("Works", "Issues", "Fails"):
        out.append(f"Registry 'results' field: {p['works']}")
    if p["uses"]:
        out.append(f"Used in {p['uses']} other parts")
    if c["verified"]:
        v = c["verified"].efficiency
        out.append(f"Verified strength value: {v.number:g} {v.unit} (legacy Registry XML)")
    if c["hint"]:
        line = f"Description says '{c['hint']['word']}' ({c['hint']['level']}), which is text, not a measurement"
        if target:
            line += " and matches your target" if c["match"] else " and does not match your target"
        out.append(line)
    elif target:
        out.append("No stated strength, so it was not matched to your target level")
    return out


def find_parts(host: str, standard: str = "RFC10", level: str | None = None, top: int = 5,
               kinds: tuple[str, ...] = KINDS, avoid_type_iis: bool = False) -> dict:
    if standard != "RFC10":
        raise PartFinderError("Only the RFC10 assembly standard is supported.")
    if level is not None and level not in LEVELS:
        raise PartFinderError(f"level must be one of {', '.join(LEVELS)}")
    if not 1 <= top <= 10:
        raise PartFinderError("top must be between 1 and 10")
    bad = [k for k in kinds if k not in KINDS]
    if bad:
        raise PartFinderError(f"Unknown part kind(s): {bad}")
    if host not in HOST_CHASSIS:
        raise PartFinderError(f"No Registry chassis tag is mapped for host {host!r}. Supported: {', '.join(sorted(HOST_CHASSIS))}")

    data, strengths = load_dataset(), load_verified_strengths()
    tag = HOST_CHASSIS[host]
    result: dict = {"host": host, "chassis_tag": tag, "standard": standard, "target_level": level, "results": {}}
    for kind in kinds:
        cands = [c for p in data["parts"] if p["kind"] == kind
                 if (c := _candidate(p, tag, level, strengths, avoid_type_iis))]
        cands.sort(key=_sort_key)
        result["results"][kind] = {
            "considered": sum(1 for p in data["parts"] if p["kind"] == kind), "passing_filters": len(cands),
            "parts": [{
                "name": c["part"]["name"], "registry_url": registry_url(c["part"]["name"]), "status": c["part"]["status"],
                "description": c["part"]["short_desc"], "chassis_evidence": c["evidence"],
                "expression": {
                    "estimate": None,
                    "verified_strength": ({"value": c["verified"].efficiency.number, "unit": c["verified"].efficiency.unit,
                                           "source": "legacy Registry XML saved 2026-09-21"} if c["verified"] else None),
                    "level_hint": c["hint"],
                    "note": ("No numeric expression estimate: the Registry holds no verified strength for this part."
                             if not c["verified"] else "Verified value is not comparable to parts measured in other units."),
                },
                "reasons": _explain(c, level),
            } for c in cands[:top]],
        }
    result["warnings"] = _warnings(result["results"])
    result["dataset"] = {"source": data["source"], "parts_in_extract": data["count"]}
    result["limitations"] = [
        "Ordering uses availability, a stated level in the description, host tag, Registry results and community use. It is not a strength prediction.",
        "Junctions between parts are not scanned for new cut sites; run the assembled sequence through the RFC10 check.",
    ]
    return result


def _warnings(results: dict) -> list[str]:
    def top_hint(kind):
        parts = results.get(kind, {}).get("parts", [])
        return parts[0]["expression"]["level_hint"] if parts else None
    p, r = top_hint("promoter"), top_hint("rbs")
    if p and r and p["level"] == "high" and r["level"] == "high":
        return ["The top promoter and top RBS are both described as strong. Very strong pairs can burden cells. "
                "This comes from description text, not measured strengths."]
    return []


def check_cds(cds: str) -> dict:
    """RFC10 status of the user's CDS (Part 1 scanner), so the finder can say whether the gene itself is compatible."""
    r = validate_rfc10(cds)
    return {"ok": r.passed, "illegal_sites": [c.enzyme for c in r.hard_failures], "warnings": [c.enzyme for c in r.warnings]}
