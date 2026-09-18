"""Amino acid -> DNA reverse translation with frequency-weighted codon
optimization.

Rather than always emitting the single most-common codon for each amino
acid (which produces unnaturally repetitive codon usage and can starve
the corresponding tRNA pool, contributing to ribosome stalling), each
codon is drawn at random with probability equal to its real usage
frequency in the target chassis. Over a full sequence this reproduces
the chassis's natural codon usage distribution.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .codon_usage import CodonUsageTable, DEFAULT_ORGANISM, load_codon_table

# Standard 20 amino acids, one-letter IUPAC codes. '*' is accepted as an
# explicit stop marker but is not part of this set since it is handled
# separately (see reverse_translate()).
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
STOP_SYMBOL = "*"


class ProteinSequenceError(ValueError):
    pass


@dataclass(frozen=True)
class CodonChoice:
    position: int  # 1-indexed amino acid position in the protein
    amino_acid: str
    codon: str
    fraction: float


@dataclass(frozen=True)
class ReverseTranslationResult:
    protein: str
    dna: str
    organism: str
    stop_codon_appended: bool
    choices: tuple[CodonChoice, ...]


def clean_protein_sequence(raw: str) -> str:
    """Strip whitespace/newlines and uppercase a pasted protein sequence."""

    return "".join(raw.split()).upper()


def validate_protein_sequence(protein: str) -> None:
    if not protein:
        raise ProteinSequenceError("Protein sequence is empty")

    body = protein[:-1] if protein.endswith(STOP_SYMBOL) else protein
    if not body:
        raise ProteinSequenceError("Protein sequence contains only a stop symbol")

    invalid = sorted(set(body) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise ProteinSequenceError(
            "Protein sequence contains characters that are not one of the "
            f"20 standard amino acids: {', '.join(invalid)!r}"
        )

    interior_stop = STOP_SYMBOL in body
    if interior_stop:
        raise ProteinSequenceError(
            f"Stop symbol {STOP_SYMBOL!r} may only appear at the end of the "
            "protein sequence"
        )


def reverse_translate(
    protein: str,
    organism: str = DEFAULT_ORGANISM,
    seed: int | None = None,
    rng: random.Random | None = None,
) -> ReverseTranslationResult:
    """Reverse-translate a protein sequence into a frequency-weighted,
    codon-optimized DNA CDS for the given chassis organism.

    A trailing stop codon is appended automatically (chosen with the same
    frequency-weighted method) if the input does not already end in '*'.

    Pass an existing ``rng`` (instead of ``seed``) to share one random
    stream with a later step, e.g. illegal-site repair -- this keeps a run
    reproducible end-to-end under a single seed.
    """

    protein = clean_protein_sequence(protein)
    validate_protein_sequence(protein)

    table = load_codon_table(organism)
    if rng is None:
        rng = random.Random(seed)

    stop_codon_appended = not protein.endswith(STOP_SYMBOL)
    residues = protein if protein.endswith(STOP_SYMBOL) else protein + STOP_SYMBOL

    codons: list[str] = []
    choices: list[CodonChoice] = []
    for position, amino_acid in enumerate(residues, start=1):
        codon, fraction = _weighted_choice(table, amino_acid, rng)
        codons.append(codon)
        choices.append(
            CodonChoice(
                position=position,
                amino_acid=amino_acid,
                codon=codon,
                fraction=fraction,
            )
        )

    dna = "".join(codons)

    return ReverseTranslationResult(
        protein=protein,
        dna=dna,
        organism=organism,
        stop_codon_appended=stop_codon_appended,
        choices=tuple(choices),
    )


def _weighted_choice(
    table: CodonUsageTable, amino_acid: str, rng: random.Random
) -> tuple[str, float]:
    options = table.codons_for(amino_acid)
    weights = [opt.fraction for opt in options]
    chosen = rng.choices(options, weights=weights, k=1)[0]
    return chosen.codon, chosen.fraction
