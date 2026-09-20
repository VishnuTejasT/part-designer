"""Report building/formatting for the multi-host optimization engine,
following report.py's ``add(...)`` banner-section text convention."""

from __future__ import annotations

from .codon_usage import load_codon_table
from .optimize import HostModeResult, Limits
from .trna_usage import TrnaDataError, load_trna_table

IN_SILICO_CAVEAT = (
    "In silico ranking cannot predict which variant performs best in vitro "
    "-- treat these as a shortlist, not a verdict."
)


def _translate(dna: str, table) -> str:
    return "".join(table.by_codon[dna[i : i + 3]]["amino_acid"] for i in range(0, len(dna) - len(dna) % 3, 3))


def build_checks(r: HostModeResult, limits: Limits) -> list[dict]:
    """The plain-English checklist rows as data: id, group, status
    (pass | review | fail), applicable, and the numbers behind it. All
    wording lives in the client's strings file -- nothing user-facing here."""

    table = load_codon_table(r.host)
    ok = r.constraint_pass_fail
    is_bacterial = table.metadata.get("domain") == "prokaryote"

    translated = _translate(r.dna, table)
    expected = r.protein if r.protein.endswith("*") else r.protein + "*"
    identity_ok = translated == expected
    # Count amino acids only; the stop codon is judged by the start/stop check.
    body_t, body_e = translated.rstrip("*"), expected.rstrip("*")
    matched = sum(1 for a, b in zip(body_t, body_e) if a == b)

    def status(passed: bool) -> str:
        return "pass" if passed else "fail"

    start_dg = r.five_prime_dG
    if start_dg is None and r.init_region_unpaired is None:
        start_status = "review"  # no start-of-gene DNA given: only an estimate is possible
    else:
        start_status = status(
            r.init_region_unpaired is not False and (start_dg is None or start_dg >= limits.init_dg_floor)
        )

    return [
        {"id": "protein_identity", "group": "must_pass", "status": status(identity_ok), "applicable": True,
         "value": {"matched": matched, "total": len(body_e)}},
        {"id": "start_stop", "group": "must_pass", "status": status(ok["no_internal_stop"]), "applicable": True,
         "value": {"ends_with_stop": translated.endswith("*"), "early_stops": translated[:-1].count("*")}},
        {"id": "cut_sites", "group": "must_pass", "status": status(ok["no_forbidden_sites"]), "applicable": True,
         "value": {"found": sum(r.forbidden_site_hits.values()), "sites": r.forbidden_site_hits}},
        {"id": "rare_codons", "group": "must_pass", "status": status(ok["zero_rare_codons"]), "applicable": True,
         "value": {"found": r.codons_below_w_threshold, "lowest": r.min_w_used, "cutoff": limits.rare_codon_w}},
        {"id": "codon_score", "group": "quality", "status": status(ok["cai_floor_met"]), "applicable": True,
         "value": {"cai": r.cai, "goal": limits.cai_floor}},
        {"id": "synthesis", "group": "quality",
         "status": status(ok["no_repeated_15mers"] and ok["no_long_homopolymers"] and ok["gc_windows_in_bounds"]),
         "applicable": True,
         "value": {"repeats": r.repeated_15mers, "gc_min": r.gc_window_min, "gc_max": r.gc_window_max}},
        {"id": "start_region", "group": "quality", "status": start_status, "applicable": True,
         "value": {"dG": start_dg, "open": r.init_region_unpaired, "limit": limits.init_dg_floor}},
        {"id": "hairpins", "group": "quality", "status": status(ok["no_overlong_stem"] and ok["no_overstable_window"]),
         "applicable": True,
         "value": {"longest_stem": r.longest_stem, "strongest": r.worst_window_mfe,
                   "stem_limit": limits.max_stem_bp - 1, "energy_limit": limits.max_window_mfe}},
        {"id": "terminators", "group": "quality", "status": status(ok["no_terminator_motifs"]),
         "applicable": is_bacterial, "value": {"found": len(r.terminator_motifs)}},
        {"id": "mirror_repeats", "group": "quality", "status": status(ok["no_inverted_repeats"]),
         "applicable": True, "value": {"found": len(r.inverted_repeats)}},
    ]


def _result_dict(r: HostModeResult, limits: Limits) -> dict:
    checks = build_checks(r, limits)
    counted = [c for c in checks if c["applicable"]]
    return {
        "checks": checks,
        "checks_total": len(counted),
        "checks_passed": sum(1 for c in counted if c["status"] == "pass"),
        "conflicts": r.conflicts,
        "host": r.host,
        "mode": r.mode,
        "protein_length": len(r.protein) - (1 if r.protein.endswith("*") else 0),
        "dna": r.dna,
        "dna_length": len(r.dna),
        "cai": r.cai,
        "cai_note": r.cai_note,
        "expression_weighted_cai": r.expression_weighted_cai,
        "expression_weighted_cai_note": r.expression_weighted_cai_note,
        "tai": r.tai,
        "tai_note": r.tai_note,
        "min_w_used": r.min_w_used,
        "codons_below_w_threshold": r.codons_below_w_threshold,
        "five_prime_dG": r.five_prime_dG,
        "five_prime_dG_note": r.five_prime_dG_note,
        "init_region_unpaired": r.init_region_unpaired,
        "worst_window_mfe": r.worst_window_mfe,
        "longest_stem": r.longest_stem,
        "terminator_motifs": r.terminator_motifs,
        "inverted_repeats": r.inverted_repeats,
        "gc_overall": r.gc_overall,
        "gc_window_min": r.gc_window_min,
        "gc_window_max": r.gc_window_max,
        "repeated_15mers": r.repeated_15mers,
        "forbidden_site_hits": r.forbidden_site_hits,
        "motif_hits": r.motif_hits,
        "pause_site_count": r.pause_site_count,
        "constraint_pass_fail": r.constraint_pass_fail,
        "overall_pass": all(r.constraint_pass_fail.values()),
        "notes": r.notes,
    }


def build_optimization_report(results: list[HostModeResult], limits: Limits | None = None) -> dict:
    limits = limits or Limits()
    hosts = sorted({r.host for r in results})

    codon_table_sources = {}
    trna_table_sources = {}
    for host in hosts:
        table = load_codon_table(host)
        codon_table_sources[host] = {
            "organism": table.metadata.get("organism"),
            "source": table.metadata.get("source"),
            "source_url": table.metadata.get("source_url"),
            "retrieved": table.metadata.get("retrieved"),
        }
        try:
            trna_table = load_trna_table(host)
            trna_table_sources[host] = {
                "organism": trna_table.metadata.get("organism"),
                "source": trna_table.metadata.get("source"),
                "source_url": trna_table.metadata.get("source_url"),
                "retrieved": trna_table.metadata.get("retrieved"),
            }
        except TrnaDataError:
            trna_table_sources[host] = None

    summary_table = [
        {
            "host": r.host,
            "mode": r.mode,
            "cai": r.cai,
            "tai": r.tai,
            "five_prime_dG": r.five_prime_dG,
            "pass": all(r.constraint_pass_fail.values()),
            "checks_passed": d["checks_passed"],
            "checks_total": d["checks_total"],
        }
        for r, d in zip(results, [_result_dict(x, limits) for x in results])
    ]

    return {
        "sequences": [_result_dict(r, limits) for r in results],
        "limits_used": {
            "cai_floor": limits.cai_floor, "rare_codon_w": limits.rare_codon_w,
            "longest_allowed_stem": limits.max_stem_bp - 1,
            "worst_window_mfe": limits.max_window_mfe, "start_region_dG": limits.init_dg_floor,
        },
        "summary_table": summary_table,
        "codon_table_sources": codon_table_sources,
        "trna_table_sources": trna_table_sources,
        "caveat": IN_SILICO_CAVEAT,
    }


def format_optimization_report_text(report: dict) -> str:
    lines: list[str] = []

    def add(text: str = "") -> None:
        lines.append(text)

    add("=" * 70)
    add("PLASMID DESIGN ASSISTANT -- Codon Optimization Engine Report")
    add("=" * 70)

    for seq in report["sequences"]:
        add()
        add("-" * 70)
        add(f"HOST: {seq['host']}  MODE: {seq['mode']}")
        add("-" * 70)
        add(f"Protein length: {seq['protein_length']} aa   DNA length: {seq['dna_length']} bp")
        add(f"DNA: {seq['dna']}")
        add()
        cai_str = f"{seq['cai']:.3f}" if seq["cai"] is not None else f"unavailable ({seq['cai_note']})"
        add(f"CAI: {cai_str}")
        ew_str = (
            f"{seq['expression_weighted_cai']:.3f}"
            if seq["expression_weighted_cai"] is not None
            else f"unavailable ({seq['expression_weighted_cai_note']})"
        )
        add(f"Expression-weighted CAI: {ew_str}")
        tai_str = f"{seq['tai']:.3f}" if seq["tai"] is not None else f"unavailable ({seq['tai_note']})"
        add(f"tAI: {tai_str}")
        add(f"Minimum w used: {seq['min_w_used']:.3f}   Codons with w<0.3: {seq['codons_below_w_threshold']}")
        dg_str = (
            f"{seq['five_prime_dG']:.2f} kcal/mol"
            if seq["five_prime_dG"] is not None
            else f"unavailable ({seq['five_prime_dG_note']})"
        )
        add(f"5' initiation dG (-20/+40): {dg_str}")
        init_str = "unavailable" if seq["init_region_unpaired"] is None else ("yes" if seq["init_region_unpaired"] else "NO")
        add(f"Initiation region (-15/+20) fully single-stranded: {init_str}")
        worst_mfe_str = f"{seq['worst_window_mfe']:.2f} kcal/mol" if seq["worst_window_mfe"] is not None else "unavailable"
        add(f"Worst 60nt window MFE: {worst_mfe_str}")
        add(f"Longest contiguous stem in any window: {seq['longest_stem']} bp")
        add(f"Terminator-like motifs: {len(seq['terminator_motifs'])}")
        add(f"Inverted repeats (>=8bp within 30nt): {len(seq['inverted_repeats'])}")
        add(
            f"GC%: overall {seq['gc_overall']*100:.1f}%  "
            f"window min {seq['gc_window_min']*100:.1f}%  "
            f"window max {seq['gc_window_max']*100:.1f}%"
        )
        add(f"Repeated 15-mers: {seq['repeated_15mers']}")
        add(f"Forbidden site hits: {seq['forbidden_site_hits'] or 'none'}")
        add(f"Cryptic/SD-like motif hits: {len(seq['motif_hits'])}")
        add(f"Pause sites placed: {seq['pause_site_count']}")
        add()
        add("Constraint pass/fail:")
        for name, passed in seq["constraint_pass_fail"].items():
            add(f"  [{'PASS' if passed else 'FAIL'}] {name}")
        add(f"Overall: {'PASS' if seq['overall_pass'] else 'FAIL'}")
        if seq["notes"]:
            add()
            add("Notes:")
            for note in seq["notes"]:
                add(f"  - {note}")

    add()
    add("=" * 70)
    add("SUMMARY")
    add("=" * 70)
    add(f"{'HOST':<24}{'MODE':<22}{'CAI':<8}{'tAI':<8}{'5dG':<10}{'PASS':<6}")
    for row in report["summary_table"]:
        cai = f"{row['cai']:.3f}" if row["cai"] is not None else "-"
        tai = f"{row['tai']:.3f}" if row["tai"] is not None else "-"
        dg = f"{row['five_prime_dG']:.2f}" if row["five_prime_dG"] is not None else "-"
        add(
            f"{row['host']:<24}{row['mode']:<22}{cai:<8}{tai:<8}{dg:<10}"
            f"{'PASS' if row['pass'] else 'FAIL':<6}"
        )

    add()
    add("Codon usage table sources:")
    for host, meta in report["codon_table_sources"].items():
        add(f"  {host}: {meta['source']} ({meta['source_url']}, retrieved {meta['retrieved']})")

    add()
    add("tRNA gene-copy-number table sources:")
    for host, meta in report["trna_table_sources"].items():
        if meta is None:
            add(f"  {host}: unavailable (no sourced tRNA GCN table)")
        else:
            add(f"  {host}: {meta['source']} ({meta['source_url']}, retrieved {meta['retrieved']})")

    add()
    add(report["caveat"])

    return "\n".join(lines)
