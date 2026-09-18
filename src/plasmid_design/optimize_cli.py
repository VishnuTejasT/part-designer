"""Command-line entry point for the multi-host codon optimization engine
(``plasmid-design-optimize``). Separate from ``cli.py``/``plasmid-design``
(Part 1 sequence converter), which is untouched."""

from __future__ import annotations

import argparse
import json
import sys

from .codon_usage import CodonUsageError, ORGANISM_TABLES
from .optimize import (
    MODES,
    OptimizationError,
    OptimizationRequest,
    StructuralRegion,
    optimize_cds,
)
from .optimize_report import build_optimization_report, format_optimization_report_text
from .reverse_translate import ProteinSequenceError


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plasmid-design-optimize",
        description=(
            "Multi-host, mode-driven codon optimization: CAI/tAI, mRNA-"
            "folding-aware initiation design, and hard sequence constraints."
        ),
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--protein", help="Amino acid sequence, one-letter codes.")
    source.add_argument(
        "--input-file",
        help="Path to a file with the raw amino acid sequence. Falls back to "
        "stdin if neither this nor --protein is given.",
    )
    parser.add_argument(
        "--hosts",
        required=True,
        help=f"Comma-separated host chassis. Choices: {', '.join(sorted(ORGANISM_TABLES))}",
    )
    parser.add_argument(
        "--mode",
        default="BALANCED",
        choices=("PRODUCTION", "FOLDING", "BALANCED", "ALL"),
        help="Optimization mode (default: BALANCED).",
    )
    parser.add_argument("--hedge", action="store_true", help="Also return a diversified variant.")
    parser.add_argument("--native-organism", default=None, help="Native source organism, for FOLDING harmonization.")
    parser.add_argument("--five-prime-utr", default="", help="5' UTR/RBS/Kozak sequence upstream of the CDS.")
    parser.add_argument("--n-tag", default="", help="N-terminal tag amino acid sequence, fused before optimization.")
    parser.add_argument("--c-tag", default="", help="C-terminal tag amino acid sequence, fused before optimization.")
    parser.add_argument(
        "--forbidden-enzymes",
        default="",
        help="Comma-separated Name:PATTERN pairs, e.g. 'EcoRV:GATATC,XhoI:CTCGAG'.",
    )
    parser.add_argument(
        "--structural-regions",
        default=None,
        help='JSON list, e.g. \'[{"start":40,"end":48,"kind":"linker"}]\'.',
    )
    parser.add_argument("--structural-regions-file", default=None, help="Path to a JSON file with the same shape as --structural-regions.")
    parser.add_argument(
        "--temperature", type=float, default=37.0,
        help="RNA folding temperature in Celsius (default: 37; use the induction temperature if known).",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--json-out", help="Optional path to also write the full structured report as JSON.")
    parser.add_argument("--quiet", action="store_true", help="Suppress the human-readable report on stdout.")
    return parser


def _read_protein_input(args: argparse.Namespace) -> str:
    if args.protein is not None:
        return args.protein
    if args.input_file is not None:
        with open(args.input_file) as fh:
            return fh.read()
    return sys.stdin.read()


def _parse_forbidden_enzymes(raw: str) -> tuple[tuple[str, str], ...]:
    if not raw.strip():
        return ()
    pairs = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid --forbidden-enzymes entry {item!r}, expected Name:PATTERN")
        name, pattern = item.split(":", 1)
        pairs.append((name.strip(), pattern.strip().upper()))
    return tuple(pairs)


def _parse_structural_regions(args: argparse.Namespace) -> tuple[StructuralRegion, ...]:
    raw = args.structural_regions
    if args.structural_regions_file:
        with open(args.structural_regions_file) as fh:
            raw = fh.read()
    if not raw:
        return ()
    data = json.loads(raw)
    return tuple(StructuralRegion(start=d["start"], end=d["end"], kind=d["kind"]) for d in data)


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    raw_protein = _read_protein_input(args)
    hosts = tuple(h.strip() for h in args.hosts.split(",") if h.strip())

    try:
        forbidden = _parse_forbidden_enzymes(args.forbidden_enzymes)
        regions = _parse_structural_regions(args)
        request = OptimizationRequest(
            protein=raw_protein,
            hosts=hosts,
            mode=args.mode,
            native_organism=args.native_organism,
            five_prime_utr=args.five_prime_utr,
            n_tag=args.n_tag,
            c_tag=args.c_tag,
            forbidden_enzymes=forbidden,
            structural_regions=regions,
            hedge=args.hedge,
            seed=args.seed,
            temperature_c=args.temperature,
        )
        results = optimize_cds(request)
    except (ProteinSequenceError, CodonUsageError, OptimizationError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    report = build_optimization_report(results)

    if not args.quiet:
        print(format_optimization_report_text(report))

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(report, fh, indent=2)

    return 0 if all(row["pass"] for row in report["summary_table"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
