import pytest
from contigbuilder.sequenceUtils import prefix, suffix, RvsComp, GetMismatches

def test_prefix_no_window():
    seq = "ACGTCATGCA"
    assert prefix(seq, len_kmer=4, n=0) == "ACGT"

def test_prefix_with_window():
    seq = "ACGTCATGCA"
    assert prefix(seq, len_kmer=4, n=2) == "GTCA"

def test_suffix_no_window():
    seq = "ACGTCATGCA"
    assert suffix(seq, len_kmer=4, n=0) == "TGCA"

def test_suffix_with_window():
    seq = "ACGTCATGCA"
    assert suffix(seq, len_kmer=4, n=2) == "CATG"

def test_reverse_complement():
    seq = "ACGTCA"
    assert RvsComp(seq) == "TGACGT"

def test_get_mismatches_zero():
    assert GetMismatches("ACGT", "ACGT") == 0

def test_get_mismatches_two():
    assert GetMismatches("ACGT", "AGGA") == 2

def test_get_mismatches_unequal_length_raises_error():
    with pytest.raises(ValueError):
        GetMismatches("ACGT", "ACG")