"""Report how many candidate parts have a valid `efficiency` value.

usage: python scripts/efficiency_coverage.py <dir-of-saved-part-XML> [candidates.csv]

Prints valid / NA-only / absent counts overall and per kind, the units seen,
and -- the number that decides whether ranking is viable -- how many of the
candidate parts have NO XML file yet (so their status is unknown, not absent).
"""
import csv
import sys

from plasmid_design.registry_params import coverage, load_xml_dir

xml_dir = sys.argv[1]
cand_csv = sys.argv[2] if len(sys.argv) > 2 else "data/registry/efficiency_probe_candidates.csv"
kinds = {r["part_name"]: r["kind"] for r in csv.DictReader(open(cand_csv))}
records = load_xml_dir(xml_dir)
c = coverage(records, kinds)
have = {r.part_name for r in records}
print(f"XML files parsed: {c['total']}   candidates in list: {len(kinds)}   candidates with no XML yet: {len(set(kinds) - have)}")
print("efficiency status:", c["status"])
for kind, s in sorted(c["by_kind"].items()):
    print(f"  {kind:11s} valid={s['valid']:4d}  na_only={s['na_only']:4d}  absent={s['absent']:4d}   units: {c['kind_units'].get(kind, {})}")
flagged = [r.part_name for r in records if r.newer_non_numeric_after_valid]
print("valid values with a newer NA row after them:", len(flagged), flagged[:5])
