"""Automatic repair of RFC10 hard-fail illegal sites.

When an illegal restriction site turns up inside a CDS, it is almost always
possible to remove it without touching the protein: most amino acids have
several synonymous codons, so re-rolling just the codon(s) that overlap the
site -- rather than the whole sequence -- can break up the illegal pattern
while leaving every translated residue identical.
"""

from __future__ import annotations

import dataclasses
import random
from dataclasses import dataclass

from .codon_usage import load_codon_table
from .reverse_translate import ReverseTranslationResult
from .rfc10 import RFC10Report, validate_rfc10

MAX_REPAIR_ITERATIONS = 200


@dataclass(frozen=True)
class CodonChange:
    position: int  # 1-indexed amino acid position that was changed
    amino_acid: str
    old_codon: str
    new_codon: str
    enzyme: str  # which illegal site this change was made to break up


@dataclass(frozen=True)
class RepairResult:
    translation: ReverseTranslationResult
    validation: RFC10Report
    changes: tuple[CodonChange, ...]
    resolved: bool
    iterations: int


def _codon_positions_overlapping(start: int, end: int) -> list[int]:
    """1-indexed amino acid positions whose codon overlaps a DNA span given
    as 1-indexed inclusive base coordinates."""

    first = (start - 1) // 3 + 1
    last = (end - 1) // 3 + 1
    return list(range(first, last + 1))


def repair_illegal_sites(
    translation: ReverseTranslationResult,
    rng: random.Random,
    max_iterations: int = MAX_REPAIR_ITERATIONS,
) -> RepairResult:
    """Resample only the codons overlapping a hard-fail RFC10 site -- never
    the amino acid sequence -- until no hard-fail sites remain or no further
    synonymous change is possible (e.g. the site falls entirely within an
    invariant codon like Met/Trp)."""

    table = load_codon_table(translation.organism)
    codons = [c.codon for c in translation.choices]
    choices = list(translation.choices)
    changes: list[CodonChange] = []

    validation = validate_rfc10("".join(codons))
    iterations = 0
    while validation.hard_failures and iterations < max_iterations:
        iterations += 1

        position_enzyme: dict[int, str] = {}
        for check in validation.hard_failures:
            for occ in check.occurrences:
                for pos in _codon_positions_overlapping(occ.start, occ.end):
                    position_enzyme.setdefault(pos, check.enzyme)

        changed_any = False
        for pos, enzyme in position_enzyme.items():
            amino_acid = choices[pos - 1].amino_acid
            options = table.codons_for(amino_acid)
            current_codon = codons[pos - 1]
            alternatives = [opt for opt in options if opt.codon != current_codon]
            if not alternatives:
                continue  # invariant amino acid (e.g. Met, Trp) -- can't change

            weights = [opt.fraction for opt in alternatives]
            if sum(weights) <= 0:
                weights = [1.0] * len(alternatives)
            chosen = rng.choices(alternatives, weights=weights, k=1)[0]

            changes.append(
                CodonChange(
                    position=pos,
                    amino_acid=amino_acid,
                    old_codon=current_codon,
                    new_codon=chosen.codon,
                    enzyme=enzyme,
                )
            )
            codons[pos - 1] = chosen.codon
            choices[pos - 1] = dataclasses.replace(
                choices[pos - 1], codon=chosen.codon, fraction=chosen.fraction
            )
            changed_any = True

        validation = validate_rfc10("".join(codons))
        if not changed_any:
            break

    repaired_translation = dataclasses.replace(
        translation, dna="".join(codons), choices=tuple(choices)
    )

    return RepairResult(
        translation=repaired_translation,
        validation=validation,
        changes=tuple(changes),
        resolved=not validation.hard_failures,
        iterations=iterations,
    )
