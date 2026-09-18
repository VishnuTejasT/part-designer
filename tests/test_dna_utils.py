from plasmid_design.dna_utils import (
    gc_content,
    gc_window_scan,
    homopolymer_runs,
    reverse_complement,
    repeated_kmers,
)


def test_reverse_complement():
    assert reverse_complement("GAATTC") == "GAATTC"
    assert reverse_complement("AAAACCC") == "GGGTTTT"


def test_gc_content():
    assert gc_content("GGCC") == 1.0
    assert gc_content("AATT") == 0.0
    assert gc_content("") == 0.0


def test_homopolymer_runs_detects_long_run_only():
    seq = "ATG" + "AAAAAAA" + "TAA"  # 7 A's, over the default max of 6
    runs = homopolymer_runs(seq)
    assert len(runs) == 1
    assert runs[0].base == "A"
    assert runs[0].length == 7


def test_homopolymer_runs_ignores_short_runs():
    seq = "ATGAAACCCTAA"  # max run length 3
    assert homopolymer_runs(seq) == []


def test_repeated_kmers_detects_duplicate_15mer():
    kmer = "ATGCATGCATGCATG"  # 15 bases
    seq = kmer + "GGG" + kmer
    reps = repeated_kmers(seq, k=15)
    assert len(reps) == 1
    assert reps[0].kmer == kmer
    assert len(reps[0].positions) == 2


def test_gc_window_scan_flags_low_and_high():
    low_seq = "A" * 60  # 0% GC
    assert gc_window_scan(low_seq, window=50)
    high_seq = "G" * 60  # 100% GC
    assert gc_window_scan(high_seq, window=50)
    mixed = "ATGC" * 15  # 50% GC, within default bounds
    assert gc_window_scan(mixed, window=50) == []
