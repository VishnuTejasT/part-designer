# Plasmid Design Assistant

A protein-to-optimized-plasmid design tool built for an iGEM Software track
contribution. It takes a team's protein sequence all the way to a
compatibility-checked, expression-optimized plasmid design.

## Part 1: Sequence Converter (implemented)

Reverse-translates a protein sequence into a codon-optimized DNA CDS for a
target chassis organism, then validates the result against the BioBrick
RFC10 assembly standard.

- **Codon optimization** is frequency-weighted: each codon is drawn at
  random with probability equal to its real usage frequency in the chassis
  organism, not always the single most common codon. This reproduces the
  organism's natural codon usage distribution instead of producing an
  unnaturally repetitive sequence that can starve a single tRNA pool.
- **Codon usage data** is real, sourced data, not a hardcoded guess: the
  default E. coli K-12 table comes from the
  [Kazusa Codon Usage Database](http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333&aa=1&style=N)
  (see `data/ecoli_k12_codon_usage.json` for the full citation and every
  codon's fraction/per-thousand/count values). The lookup is organized so a
  new chassis organism can be added by dropping in a new JSON file and
  registering it in `ORGANISM_TABLES` in `codon_usage.py` -- nothing else in
  the codebase is organism-specific.
- **RFC10 validation** checks both strands (forward and reverse complement)
  for the four illegal restriction sites (EcoRI, XbaI, SpeI, PstI -- hard
  fail) and flags NotI as a warning, per the confirmed BioBrick RFC10 rule
  set. Every conflict is reported with its exact position and strand, not
  just pass/fail. Golden Gate / BsaI-BsmBI overhang validation is
  intentionally out of scope for now.

### Usage

```bash
pip install -e .

plasmid-design --protein "MKTAYIAKQRQISFVKSHFSRQ" --seed 42
```

Or read the sequence from a file or stdin:

```bash
plasmid-design --input-file my_protein.txt
echo "MKTAYIAKQRQISFVKSHFSRQ" | plasmid-design
```

Write the full structured report to JSON as well as stdout:

```bash
plasmid-design --protein "MKTAYIAKQRQISFVKSHFSRQ" --json-out report.json
```

`--seed` makes codon selection reproducible; omit it for a fresh random draw
each run. Exit code is `0` if the sequence passes RFC10 (hard-fail sites
only), `1` if it fails, `2` on invalid input.

### Tests

```bash
pip install -e ".[dev]"
pytest
```

## Part 2: Part Finder (not yet started)

Will recommend compatible flanking parts (promoter, RBS, terminator) from
the real iGEM Part Registry for the validated CDS from Part 1.

## Project layout

```
data/                          Sourced codon usage tables (JSON), one per chassis organism
src/plasmid_design/
    codon_usage.py              Loads/validates a chassis's codon usage table
    reverse_translate.py        Protein -> frequency-weighted codon-optimized DNA
    rfc10.py                    BioBrick RFC10 illegal-site scanning (both strands)
    report.py                   Structured + human-readable report assembly
    cli.py                      Command-line entry point
tests/                          pytest unit tests for every module above
```
