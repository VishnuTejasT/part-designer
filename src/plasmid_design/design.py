"""Top-level pipeline: protein -> codon-optimized CDS -> RFC10 validation,
with automatic illegal-site repair by default."""

from __future__ import annotations

import random
from dataclasses import dataclass

from .codon_usage import DEFAULT_ORGANISM
from .repair import MAX_REPAIR_ITERATIONS, RepairResult, repair_illegal_sites
from .reverse_translate import ReverseTranslationResult, reverse_translate
from .rfc10 import RFC10Report, validate_rfc10


@dataclass(frozen=True)
class DesignResult:
    translation: ReverseTranslationResult
    validation: RFC10Report
    repair: RepairResult | None  # None if auto_fix was off or nothing needed fixing


def design_cds(
    protein: str,
    organism: str = DEFAULT_ORGANISM,
    seed: int | None = None,
    auto_fix: bool = True,
    max_repair_iterations: int = MAX_REPAIR_ITERATIONS,
) -> DesignResult:
    """Reverse-translate ``protein`` and validate it against RFC10. If
    ``auto_fix`` is set (the default) and hard-fail illegal sites are found,
    attempt to remove them by resampling only the overlapping codons -- the
    protein sequence never changes.

    A single random stream (seeded from ``seed``) drives both the initial
    codon selection and any repair, so the whole run is reproducible.
    """

    rng = random.Random(seed)
    translation = reverse_translate(protein, organism=organism, rng=rng)
    validation = validate_rfc10(translation.dna)

    repair: RepairResult | None = None
    if auto_fix and validation.hard_failures:
        repair = repair_illegal_sites(
            translation, rng, max_iterations=max_repair_iterations
        )
        translation = repair.translation
        validation = repair.validation

    return DesignResult(translation=translation, validation=validation, repair=repair)
