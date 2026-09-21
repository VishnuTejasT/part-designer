"""Part 3: assemble backbone + promoter + RBS + CDS + terminator (BioBrick RFC10)
and check the finished plasmid, not just the pieces.

Reuses Part 1's scanner (rfc10.validate_rfc10 / forbidden_sites.scan_forbidden_sites).
Part 1 had no scar logic, so the two scars are defined here.

Every threshold below is either cited or labeled as a project setting:
  * Scars: Registry pages Help:Standards/Assembly/RFC10 and Help:Assembly/Scars.
  * Synthesis limits: Twist Bioscience gene-synthesis FAQ.
  * Plasmid size: a study reports lower transformation efficiency as pUC-derived
    plasmids grow (2.6-16.1 kb tested). It gives a trend, not a cutoff, so size is
    reported with that evidence and NOT turned into an invented pass/fail line.
"""

from __future__ import annotations

import json
from functools import lru_cache

from .codon_usage import DATA_DIR
from .dna_utils import gc_content
from .forbidden_sites import scan_forbidden_sites
from .part_finder import load_dataset
from .rfc10 import validate_rfc10

SCAR_STANDARD = "TACTAGAG"   # 8 bp, between standard parts (frame-shifting in coding regions)
SCAR_RBS_CDS = "TACTAG"      # 6 bp, between an RBS and a coding sequence
HARD_ENZYMES = ("EcoRI", "XbaI", "SpeI", "PstI")

TWIST_WINDOW_GC = (10.0, 90.0)   # a 50 bp window outside this is 'high complexity' at Twist
TWIST_CLONAL_MAX_BP = 7000       # clonal gene length range 300-7,000 bp
PROJECT_GC_RANGE = (30.0, 70.0)  # this project's own design range (Part 1); NOT a vendor limit
SIZE_STUDY_MAX_BP = 16100        # largest size in the cited pUC-derived size study

SOURCES = {
    "scars": ["https://parts.igem.org/Help:Standards/Assembly/RFC10", "https://parts.igem.org/Help:Assembly/Scars"],
    "synthesis": "https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow",
    "plasmid_size": "The Effect of Increasing Plasmid Size on Transformation Efficiency in Escherichia coli "
                    "(pUC-derived plasmids, 2.6-16.1 kb): https://www.researchgate.net/publication/267219864",
}


class VectorizerError(ValueError):
    pass


@lru_cache(maxsize=1)
def load_backbones() -> dict:
    with open(DATA_DIR / "registry" / "backbones.json") as fh:
        doc = json.load(fh)
    return {b["name"]: {**b, "source": doc["source"]} for b in doc["backbones"]}


def _part(name: str, kind: str) -> dict:
    for p in load_dataset()["parts"]:
        if p["name"] == name:
            if p["kind"] != kind:
                raise VectorizerError(f"{name} is a {p['kind']}, not a {kind}.")
            return p
    raise VectorizerError(f"{name} is not in the Registry extract.")


def _hits(seq: str) -> set:
    """(enzyme, start, end) of every RFC10 hard-fail site, scanning the circular sequence."""
    circ = seq + seq[:5]
    return {(c.enzyme, o.start, o.end) for c in validate_rfc10(circ).hard_failures for o in c.occurrences}


def assemble(promoter: str, rbs: str, cds: str, terminator: str, backbone: str = "pSB1C3",
             avoid_type_iis: bool = False) -> dict:
    bbs = load_backbones()
    if backbone not in bbs:
        raise VectorizerError(f"Unknown backbone {backbone!r}. Available: {', '.join(sorted(bbs))}")
    cds = "".join(cds.split()).upper()
    if not cds or set(cds) - set("ACGT"):
        raise VectorizerError("cds must be DNA (A, C, G, T only).")
    if len(cds) % 3 or not cds.startswith("ATG") or cds[-3:] not in ("TAA", "TAG", "TGA"):
        raise VectorizerError("cds must start with ATG, end with a stop codon and be a whole number of codons.")

    B = bbs[backbone]["sequence"]
    pieces = [("promoter", promoter, _part(promoter, "promoter")["sequence"]), ("scar", None, SCAR_STANDARD),
              ("rbs", rbs, _part(rbs, "rbs")["sequence"]), ("scar", None, SCAR_RBS_CDS),
              ("cds", "your CDS", cds), ("scar", None, SCAR_STANDARD),
              ("terminator", terminator, _part(terminator, "terminator")["sequence"])]
    insert, layout, pos = "", [], len(B)
    for kind, label, seq in pieces:
        layout.append({"kind": kind, "name": label, "start": pos + 1, "end": pos + len(seq), "length": len(seq)})
        insert += seq; pos += len(seq)
    full = B + insert          # circular: the end of the insert joins the start of the backbone

    baseline = {h for h in _hits(B) if h[2] <= len(B)}       # the flank sites every BioBrick backbone carries
    new = sorted(h for h in _hits(full) if h not in baseline)
    seams = [e for x in layout for e in (x["start"], x["end"]) if x["kind"] == "scar"] + [len(B), len(B) + 1, len(full)]
    violations = [{"enzyme": e, "start": s, "end": t,
                   "at_junction": any(s - 8 <= m <= t + 8 for m in seams) or t > len(full)} for e, s, t in new]
    type_iis = sorted(n for n in scan_forbidden_sites(insert) if n in ("BsaI", "BsmBI", "BbsI", "SapI")) if avoid_type_iis else []

    win = [gc_content(insert[i:i + 50]) * 100 for i in range(max(1, len(insert) - 49))]
    lo, hi = min(win), max(win)
    gc = gc_content(insert) * 100
    total = len(full)
    return {
        "backbone": {"name": backbone, "length": len(B), "description": bbs[backbone]["description"],
                     "registry_id": bbs[backbone]["registry_id"], "source": bbs[backbone]["source"],
                     "expected_flank_sites": sorted({h[0] for h in baseline} | {"NotI"})},
        "layout": layout, "scars": {"standard": SCAR_STANDARD, "rbs_to_cds": SCAR_RBS_CDS},
        "insert": insert, "sequence": full, "fasta": f">construct_{backbone}_insert\n" + "\n".join(full[i:i + 60] for i in range(0, total, 60)) + "\n",
        "junction_check": {
            "ok": not violations and not type_iis, "violations": violations, "type_iis_sites": type_iis,
            "note": "Checks the whole circular plasmid and ignores only the sites the empty backbone already has."},
        "size": {"total_bp": total, "insert_bp": len(insert), "backbone_bp": len(B),
                 "within_studied_range": total <= SIZE_STUDY_MAX_BP,
                 "evidence": ("Transformation efficiency falls as pUC-derived plasmids get larger (studied 2.6-16.1 kb); "
                              "no cutoff is published there, so this is a trend, not a pass/fail line."
                              + ("" if total <= SIZE_STUDY_MAX_BP else " This construct is larger than the studied range.")),
                 "source": SOURCES["plasmid_size"]},
        "gc": {"insert_percent": round(gc, 1), "window_min": round(lo, 1), "window_max": round(hi, 1),
               "twist_high_complexity": lo < TWIST_WINDOW_GC[0] or hi > TWIST_WINDOW_GC[1],
               "outside_project_range": not PROJECT_GC_RANGE[0] <= gc <= PROJECT_GC_RANGE[1],
               "twist_note": "Twist calls a sequence high-complexity if a 50 bp window is below 10% or above 90% GC.",
               "project_note": "30-70% is this project's design range, not a vendor limit.", "source": SOURCES["synthesis"]},
        "synthesis": {"insert_bp": len(insert), "beyond_clonal_gene_limit": len(insert) > TWIST_CLONAL_MAX_BP,
                      "note": "Twist clonal genes are 300-7,000 bp.", "source": SOURCES["synthesis"]},
        "sources": SOURCES,
    }
