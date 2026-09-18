"""Loads chassis-specific codon usage frequency tables.

Each organism's table is a JSON file under ``data/`` containing real,
sourced codon usage bias data (see each file's ``source``/``source_url``
fields). Adding a new chassis means dropping in a new JSON file and
registering it in ``ORGANISM_TABLES`` below -- nothing else in this
package is organism-specific.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data"

# Maps an organism key (used throughout the CLI/API) to its codon usage
# data file under DATA_DIR. E. coli K-12 is the default/first chassis;
# add more organisms here as their frequency tables are sourced.
ORGANISM_TABLES: dict[str, str] = {
    "e_coli_k12": "ecoli_k12_codon_usage.json",
    "e_coli_bl21_de3": "e_coli_bl21_de3_codon_usage.json",
    "b_subtilis_168": "b_subtilis_168_codon_usage.json",
    "s_cerevisiae": "s_cerevisiae_codon_usage.json",
    "p_pastoris": "p_pastoris_codon_usage.json",
    "c_glutamicum": "c_glutamicum_codon_usage.json",
    "p_putida": "p_putida_codon_usage.json",
    "synechocystis_pcc6803": "synechocystis_pcc6803_codon_usage.json",
    "c_reinhardtii": "c_reinhardtii_codon_usage.json",
    "a_tumefaciens": "a_tumefaciens_codon_usage.json",
    "l_lactis": "l_lactis_codon_usage.json",
    "b_megaterium": "b_megaterium_codon_usage.json",
    "human": "human_codon_usage.json",
}

DEFAULT_ORGANISM = "e_coli_k12"

_FRACTION_SUM_TOLERANCE = 0.03


@dataclass(frozen=True)
class CodonOption:
    codon: str
    fraction: float
    per_thousand: float


@dataclass(frozen=True)
class CodonUsageTable:
    organism: str
    metadata: dict
    by_codon: dict[str, dict]
    by_amino_acid: dict[str, tuple[CodonOption, ...]]

    def codons_for(self, amino_acid: str) -> tuple[CodonOption, ...]:
        try:
            return self.by_amino_acid[amino_acid]
        except KeyError as exc:
            raise KeyError(
                f"No codons found for amino acid symbol {amino_acid!r} "
                f"in {self.organism} codon usage table"
            ) from exc


class CodonUsageError(ValueError):
    pass


@lru_cache(maxsize=None)
def load_codon_table(organism: str = DEFAULT_ORGANISM) -> CodonUsageTable:
    """Load (and cache) the codon usage table for a chassis organism."""

    if organism not in ORGANISM_TABLES:
        available = ", ".join(sorted(ORGANISM_TABLES))
        raise CodonUsageError(
            f"Unknown organism {organism!r}. Available chassis: {available}"
        )

    path = DATA_DIR / ORGANISM_TABLES[organism]
    if not path.exists():
        raise CodonUsageError(f"Codon usage data file missing: {path}")

    with path.open() as fh:
        raw = json.load(fh)

    by_codon: dict[str, dict] = raw["codons"]

    by_amino_acid: dict[str, list[CodonOption]] = {}
    for codon, info in by_codon.items():
        aa = info["amino_acid"]
        by_amino_acid.setdefault(aa, []).append(
            CodonOption(
                codon=codon,
                fraction=float(info["fraction"]),
                per_thousand=float(info["per_thousand"]),
            )
        )

    _validate_fractions(organism, by_amino_acid)

    frozen_by_aa = {aa: tuple(opts) for aa, opts in by_amino_acid.items()}

    metadata = {k: v for k, v in raw.items() if k != "codons"}

    return CodonUsageTable(
        organism=organism,
        metadata=metadata,
        by_codon=by_codon,
        by_amino_acid=frozen_by_aa,
    )


def _validate_fractions(
    organism: str, by_amino_acid: dict[str, list[CodonOption]]
) -> None:
    """Sanity-check that each amino acid's synonymous-codon fractions
    sum to ~1.0, catching a malformed or mistranscribed data file early."""

    for aa, options in by_amino_acid.items():
        total = sum(opt.fraction for opt in options)
        if abs(total - 1.0) > _FRACTION_SUM_TOLERANCE:
            raise CodonUsageError(
                f"Codon usage table for {organism!r} is malformed: "
                f"fractions for amino acid {aa!r} sum to {total:.3f}, "
                "expected ~1.0"
            )
