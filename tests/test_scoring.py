import pytest
from contigbuilder.scoring import bestseq

def make_hit(read_id, seq, window=0, strand="fwd"):
    return {"readId": read_id,"seq": seq,"strand": strand,"window": window}

def test_bestseq_right_window_zero():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("read1", "CCCCGGGGTTTT", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is not None
    assert best["hit"]["readId"] == "read1"
    assert best["new_contig"] == "AAAACCCCGGGGTTTT"
    assert best["overlap_len"] == 8
    assert best["mismatches"] == 0
    assert best["mismatch_rate"] == 0
    assert best["lengthAdded"] == 4

def test_bestseq_right_window_three():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("read2", "AAACCCCGGGGTTTT", window=3)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is not None
    assert best["hit"]["readId"] == "read2"
    assert best["new_contig"] == "AAAACCCCGGGGTTTT"
    assert best["overlap_len"] == 11
    assert best["mismatches"] == 0
    assert best["lengthAdded"] == 4

def test_bestseq_left_window_zero():
    contig = "CCCCGGGGTTTT"
    hits = [make_hit("read3", "AAAACCCCGGGG", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="left",len_kmer=8,mRate=0.1)

    assert best is not None
    assert best["hit"]["readId"] == "read3"
    assert best["new_contig"] == "AAAACCCCGGGGTTTT"
    assert best["overlap_len"] == 8
    assert best["mismatches"] == 0
    assert best["lengthAdded"] == 4

def test_bestseq_left_window_three():
    contig = "AAACCCCGGGGTTTT"
    hits = [make_hit("read4", "TTTTAAACCCCGGGG", window=3)]

    best = bestseq(contig=contig,hits=hits,direction="left",len_kmer=8,mRate=0.1)

    assert best is not None
    assert best["hit"]["readId"] == "read4"
    assert best["new_contig"] == "TTTTAAACCCCGGGGTTTT"
    assert best["overlap_len"] == 11
    assert best["mismatches"] == 0
    assert best["lengthAdded"] == 4

def test_bestseq_rejects_poor_candidate():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("bad_read", "TTTTGGGGAAAA", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is None

def test_bestseq_allows_one_mismatch_with_relaxed_threshold():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("one_mismatch", "CCCTGGGGTTTT", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.2)

    assert best is not None
    assert best["hit"]["readId"] == "one_mismatch"
    assert best["mismatches"] == 1
    assert best["mismatch_rate"] == 1 / 8
    assert best["new_contig"] == "AAAACCCCGGGGTTTT"

def test_bestseq_rejects_one_mismatch_with_strict_threshold():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("one_mismatch", "CCCTGGGGTTTT", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is None

def test_bestseq_chooses_highest_scoring_candidate():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("short_clean", "CCCCGGGGTT", window=0),make_hit("long_clean", "CCCCGGGGTTTTTT", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is not None
    assert best["hit"]["readId"] == "long_clean"
    assert best["new_contig"] == "AAAACCCCGGGGTTTTTT"

def test_bestseq_rejects_zero_extension():
    contig = "AAAACCCCGGGG"
    hits = [make_hit("contained_read", "CCCCGGGG", window=0)]

    best = bestseq(contig=contig,hits=hits,direction="right",len_kmer=8,mRate=0.1)

    assert best is None