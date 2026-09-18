"""Command-line entry point for Part 1: protein -> optimized, validated CDS."""

from __future__ import annotations

import argparse
import json
import sys

from .codon_usage import CodonUsageError, DEFAULT_ORGANISM, ORGANISM_TABLES
from .design import design_cds
from .report import build_report, format_report_text
from .reverse_translate import ProteinSequenceError


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plasmid-design",
        description=(
            "Reverse-translate a protein sequence into a frequency-weighted, "
            "codon-optimized DNA CDS and validate it against BioBrick RFC10."
        ),
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--protein", help="Amino acid sequence, one-letter codes (raw text)."
    )
    source.add_argument(
        "--input-file",
        help="Path to a file containing the raw amino acid sequence. "
        "If omitted and --protein is not given, the sequence is read from stdin.",
    )
    parser.add_argument(
        "--organism",
        default=DEFAULT_ORGANISM,
        choices=sorted(ORGANISM_TABLES),
        help=f"Target chassis organism (default: {DEFAULT_ORGANISM}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for codon selection, for reproducible output.",
    )
    parser.add_argument(
        "--no-auto-fix",
        action="store_true",
        help="Report illegal RFC10 sites without attempting to remove them "
        "by resampling the overlapping codon(s).",
    )
    parser.add_argument(
        "--json-out",
        help="Optional path to also write the full structured report as JSON.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the human-readable report on stdout (useful with --json-out).",
    )
    return parser


def _read_protein_input(args: argparse.Namespace) -> str:
    if args.protein is not None:
        return args.protein
    if args.input_file is not None:
        with open(args.input_file) as fh:
            return fh.read()
    return sys.stdin.read()


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    raw_protein = _read_protein_input(args)

    try:
        result = design_cds(
            raw_protein,
            organism=args.organism,
            seed=args.seed,
            auto_fix=not args.no_auto_fix,
        )
    except (ProteinSequenceError, CodonUsageError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    report = build_report(result.translation, result.validation, result.repair)

    if not args.quiet:
        print(format_report_text(report))

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(report, fh, indent=2)

    return 0 if result.validation.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
