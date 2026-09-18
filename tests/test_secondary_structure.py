import pytest

pytest.importorskip("RNA", reason="ViennaRNA not installed in this environment")

from plasmid_design.secondary_structure import (
    inverted_repeats,
    longest_stem_length,
    sliding_window_fold,
    terminator_like_motifs,
    worst_window_metrics,
)


def test_longest_stem_length_of_simple_hairpin():
    # "((((....))))" -> a single contiguous 4bp stem.
    assert longest_stem_length("((((....))))") == 4


def test_longest_stem_length_ignores_bulge_interruptions():
    # A stem broken by an interior loop should not count as one long run.
    structure = "((.((...)).))"
    assert longest_stem_length(structure) < 4


def test_sliding_window_fold_covers_short_sequence():
    windows = sliding_window_fold("ATGAAACTGGTCTAA", window=60, step=6)
    assert windows is not None
    assert len(windows) == 1
    assert windows[0].start == 1


def test_worst_window_metrics_of_stable_hairpin():
    seq = "ATG" + "GCGCGCGC" + "AAAA" + "GCGCGCGC" + "TAA"
    windows = sliding_window_fold(seq, window=40, step=4)
    worst_mfe, longest_stem = worst_window_metrics(windows)
    assert worst_mfe < 0
    assert longest_stem >= 8


def test_terminator_like_motif_detected():
    # GC-rich self-complementary stem immediately followed by a poly-T run.
    seq = "ATG" + "GCGCGCGC" + "AAAA" + "GCGCGCGC" + "TTTTT" + "AAATAA"
    windows = sliding_window_fold(seq, window=40, step=4)
    hits = terminator_like_motifs(seq, windows)
    assert len(hits) >= 1
    assert hits[0].stem_length >= 6
    assert hits[0].t_run_length >= 4


def test_no_terminator_motif_without_poly_t():
    seq = "ATG" + "GCGCGCGC" + "AAAA" + "GCGCGCGC" + "AAAAA" + "AAATAA"
    windows = sliding_window_fold(seq, window=40, step=4)
    assert terminator_like_motifs(seq, windows) == []


def test_inverted_repeat_detected_within_span():
    arm = "ACGTACGT"
    seq = "ATG" + arm + "A" * 10 + arm + "TAA"  # self-revcomp arm
    hits = inverted_repeats(seq, min_len=8, max_span=30)
    assert len(hits) >= 1


def test_no_inverted_repeat_beyond_span():
    arm = "ACGTACGT"
    seq = "ATG" + arm + "A" * 40 + arm + "TAA"  # arms too far apart
    hits = inverted_repeats(seq, min_len=8, max_span=30)
    assert hits == []


def test_no_inverted_repeat_in_random_looking_sequence():
    assert inverted_repeats("TTAGTTGTGCCGCAGCGAAGTAGTGCTTGA", min_len=8, max_span=30) == []
