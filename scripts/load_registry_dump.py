"""Load the iGEM Registry MySQL dump (parts + parts_seq_features) into SQLite.

Reads only local files; makes no network requests. Records the source file's
name, size and SHA-256 in a `provenance` table so every downstream part can be
traced to the exact dump it came from.

usage: python scripts/load_registry_dump.py "<dump.sql>" <out.sqlite>
"""
import hashlib
import re
import sqlite3
import sys
from pathlib import Path

PARTS_COLUMNS = (
    "part_id ok part_name short_desc description part_type author owning_group_id status dominant informational "
    "discontinued part_status sample_status p_status_cache s_status_cache creation_date m_datetime m_user_id uses "
    "doc_size works favorite specified_u_list deep_u_list deep_count ps_string scars default_scars owner_id "
    "group_u_list has_barcode notes source nickname categories sequence sequence_sha1 sequence_update seq_edit_cache "
    "review_result review_count review_total flag sequence_length temp_1 temp_2 temp_3 temp4 rating"
).split()
FEATURE_COLUMNS = "feature_id feature_type start_pos end_pos label part_id mark old reverse".split()
KEEP = ("part_id ok part_name short_desc description part_type status part_status sample_status discontinued "
        "creation_date m_datetime uses works categories notes source sequence sequence_length rating "
        "review_result review_count review_total").split()

TOKEN = re.compile(
    r"""(?:_binary\s+)?'((?:[^'\\]|\\.|'')*)'      # quoted string (optionally _binary)
       |(NULL)
       |(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)
       |(0x[0-9A-Fa-f]+)""",
    re.X | re.S,
)
ESC = {"n": "\n", "r": "\r", "t": "\t", "0": "\0", "b": "\b", "Z": "\x1a", "\\": "\\", "'": "'", '"': '"'}
UNESC = re.compile(r"\\(.)|''", re.S)


def unescape(s: str) -> str:
    if "\\" not in s and "''" not in s:
        return s
    return UNESC.sub(lambda m: "'" if m.group(0) == "''" else ESC.get(m.group(1), m.group(1)), s)


def rows(values_text: str, ncols: int):
    vals, out = [], []
    pos = 0
    for m in TOKEN.finditer(values_text):
        gap = values_text[pos:m.start()]
        if gap.strip(" ,()\n;") != "":
            raise ValueError(f"unparsed text between values: {gap[:60]!r}")
        pos = m.end()
        if m.group(1) is not None:
            vals.append(unescape(m.group(1)))
        elif m.group(2):
            vals.append(None)
        elif m.group(3):
            t = m.group(3)
            vals.append(float(t) if any(c in t for c in ".eE") else int(t))
        else:
            vals.append(m.group(4))
        if len(vals) == ncols:
            out.append(tuple(vals)); vals = []
    if vals:
        raise ValueError(f"leftover {len(vals)} values; column count mismatch")
    return out


def main(dump: str, out: str):
    dump_path = Path(dump)
    sha = hashlib.sha256()
    with dump_path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
    db = sqlite3.connect(out)
    db.executescript(
        "DROP TABLE IF EXISTS parts; DROP TABLE IF EXISTS features; DROP TABLE IF EXISTS provenance;"
        f"CREATE TABLE parts ({', '.join(KEEP)});"
        f"CREATE TABLE features ({', '.join(FEATURE_COLUMNS)});"
        "CREATE TABLE provenance (key TEXT, value TEXT);"
    )
    idx = [PARTS_COLUMNS.index(c) for c in KEEP]
    n_parts = n_feat = 0
    with dump_path.open("r", encoding="latin-1", newline="") as fh:
        for line in fh:
            if line.startswith("INSERT INTO `parts` VALUES "):
                body = line[len("INSERT INTO `parts` VALUES "):]
                for r in rows(body, len(PARTS_COLUMNS)):
                    db.execute(f"INSERT INTO parts VALUES ({','.join('?' * len(KEEP))})", [r[i] for i in idx]); n_parts += 1
            elif line.startswith("INSERT INTO `parts_seq_features` VALUES "):
                body = line[len("INSERT INTO `parts_seq_features` VALUES "):]
                for r in rows(body, len(FEATURE_COLUMNS)):
                    db.execute(f"INSERT INTO features VALUES ({','.join('?' * len(FEATURE_COLUMNS))})", r); n_feat += 1
    for k, v in (("source_file", dump_path.name), ("source_bytes", str(dump_path.stat().st_size)),
                 ("source_sha256", sha.hexdigest()), ("parts_rows", str(n_parts)), ("feature_rows", str(n_feat)),
                 ("dump_header", "MySQL dump, database parts_bbdb (iGEM Registry point-in-time dump)")):
        db.execute("INSERT INTO provenance VALUES (?,?)", (k, v))
    db.execute("CREATE INDEX parts_name ON parts(part_name)"); db.execute("CREATE INDEX feat_part ON features(part_id)")
    db.commit(); db.close()
    print(f"parts={n_parts} features={n_feat} sha256={sha.hexdigest()}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
