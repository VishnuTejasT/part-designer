"""Build data/registry/parts_dataset.json.gz from the SQLite made by load_registry_dump.py.

Keeps non-deleted promoters, RBS parts and terminators, with the fields Part
Finder needs, plus the dump's provenance (file name, size, SHA-256). Nothing is
added from outside the dump.

usage: python scripts/build_parts_dataset.py <registry.sqlite> [out.json.gz]
"""
import gzip
import json
import re
import sqlite3
import sys

db = sqlite3.connect(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else "data/registry/parts_dataset.json.gz"
prov = dict(db.execute("select key, value from provenance"))
KIND = "case when part_type='RBS' then 'rbs' when part_type='Terminator' then 'terminator' else 'promoter' end"
rows = db.execute(f"""
  select part_name, {KIND}, status, categories, works, uses, sequence, short_desc, part_id
  from parts
  where part_name like 'BBa_%' and status != 'Deleted' and sequence != ''
    and (part_type in ('RBS','Terminator') or (part_type='Regulatory' and categories like '%//promoter%'))
  order by part_name""").fetchall()
parts = [{
    "name": n, "kind": k, "status": s, "chassis": sorted(set(re.findall(r"//chassis/[a-z0-9_/]*[a-z0-9_]", c or ""))),
    "works": w or "", "uses": u if isinstance(u, int) and u >= 0 else None, "sequence": q.strip().upper(),
    "short_desc": (d or "")[:100], "registry_id": pid,
} for n, k, s, c, w, u, q, d, pid in rows]
doc = {"schema": 1, "source": {"kind": "iGEM Registry SQL dump (parts_bbdb)", "file": prov["source_file"],
       "bytes": int(prov["source_bytes"]), "sha256": prov["source_sha256"], "parts_rows_in_dump": int(prov["parts_rows"])},
       "count": len(parts), "parts": parts}
with gzip.open(out, "wt", encoding="utf-8", compresslevel=9) as f:
    json.dump(doc, f, separators=(",", ":"))
print(out, len(parts), "parts")
