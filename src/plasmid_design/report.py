"""Assembles and formats the Part 1 output report: DNA sequence, codon
optimization detail, and RFC10 validation results with exact conflict
positions."""

from __future__ import annotations

from .codon_usage import load_codon_table
from .repair import RepairResult
from .reverse_translate import ReverseTranslationResult
from .rfc10 import RFC10Report


def build_report(
    translation: ReverseTranslationResult,
    validation: RFC10Report,
    repair: RepairResult | None = None,
) -> dict:
    table = load_codon_table(translation.organism)

    return {
        "organism": translation.organism,
        "codon_usage_source": {
            "source": table.metadata.get("source"),
            "source_url": table.metadata.get("source_url"),
            "retrieved": table.metadata.get("retrieved"),
        },
        "protein": {
            "sequence": translation.protein,
            "length": len(translation.protein),
            "stop_codon_appended": translation.stop_codon_appended,
        },
        "dna": {
            "sequence": translation.dna,
            "length": len(translation.dna),
        },
        "codon_choices": [
            {
                "position": c.position,
                "amino_acid": c.amino_acid,
                "codon": c.codon,
                "usage_fraction": c.fraction,
            }
            for c in translation.choices
        ],
        "validation": {
            "overall_pass": validation.passed,
            "checks": [
                {
                    "enzyme": check.enzyme,
                    "pattern": check.pattern,
                    "severity": check.severity,
                    "passed": check.passed,
                    "occurrences": [
                        {
                            "strand": occ.strand,
                            "start": occ.start,
                            "end": occ.end,
                            "matched_sequence": occ.matched_sequence,
                        }
                        for occ in check.occurrences
                    ],
                }
                for check in validation.checks
            ],
        },
        "repair": {
            "attempted": repair is not None,
            "resolved": repair.resolved if repair is not None else None,
            "iterations": repair.iterations if repair is not None else 0,
            "changes": [
                {
                    "position": c.position,
                    "amino_acid": c.amino_acid,
                    "old_codon": c.old_codon,
                    "new_codon": c.new_codon,
                    "enzyme": c.enzyme,
                }
                for c in (repair.changes if repair is not None else ())
            ],
        },
    }


def format_report_text(report: dict) -> str:
    lines: list[str] = []
    add = lines.append

    add("=" * 70)
    add("PLASMID DESIGN ASSISTANT -- Part 1: Sequence Converter Report")
    add("=" * 70)

    add("")
    add(f"Chassis organism: {report['organism']}")
    src = report["codon_usage_source"]
    add(f"Codon usage source: {src['source']}")
    add(f"  {src['source_url']}")

    add("")
    add(f"Protein sequence ({report['protein']['length']} aa):")
    add(f"  {report['protein']['sequence']}")
    if report["protein"]["stop_codon_appended"]:
        add("  (no stop codon was present in the input; one was appended)")

    add("")
    add(f"Optimized DNA CDS ({report['dna']['length']} bp):")
    add(f"  {report['dna']['sequence']}")

    repair = report["repair"]
    if repair["attempted"]:
        add("")
        add("-" * 70)
        add("AUTOMATIC ILLEGAL-SITE REPAIR")
        add("-" * 70)
        n = len(repair["changes"])
        if repair["resolved"]:
            add(
                f"Resolved: changed {n} codon(s) over {repair['iterations']} "
                "pass(es); the protein sequence is unchanged."
            )
        else:
            add(
                f"NOT fully resolved after {repair['iterations']} pass(es) "
                f"({n} codon(s) changed). One or more illegal sites could not "
                "be removed by a synonymous codon change alone -- see below."
            )
        for c in repair["changes"]:
            add(
                f"  position {c['position']} ({c['amino_acid']}): "
                f"{c['old_codon']} -> {c['new_codon']}  "
                f"[breaking up a {c['enzyme']} site]"
            )

    add("")
    add("-" * 70)
    add("RFC10 COMPATIBILITY VALIDATION")
    add("-" * 70)

    validation = report["validation"]
    for check in validation["checks"]:
        label = "HARD FAIL" if check["severity"] == "fail" else "WARNING"
        status = "PASS" if check["passed"] else label
        add(
            f"[{status:>9}] {check['enzyme']:<6} "
            f"site {check['pattern']} ({check['severity']})"
        )
        for occ in check["occurrences"]:
            add(
                f"             found on strand {occ['strand']} at "
                f"position {occ['start']}-{occ['end']} "
                f"(sequence: {occ['matched_sequence']})"
            )

    add("")
    overall = "PASS" if validation["overall_pass"] else "FAIL"
    add(f"Overall RFC10 result: {overall}")
    if not validation["overall_pass"]:
        add(
            "  One or more illegal restriction sites were found in the CDS. "
            "This sequence is NOT safe to submit for BioBrick RFC10 assembly "
            "until the flagged codon(s) are re-optimized to remove them."
        )

    add("=" * 70)
    return "\n".join(lines)
