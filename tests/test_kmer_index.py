from contigbuilder.kmer_index import (createCandidate,pref_coords,suf_coords,get_NKmers_PrefSuftables)

def test_create_candidate_stores_expected_fields():
    headers = {"AAAACCCCGGGG": "read1"}

    candidate = createCandidate(read="AAAACCCCGGGG",seq="AAAACCCCGGGG",kmer="AAAA",win=0,headers=headers,strand="fwd",k_start=0,k_end=4)

    assert candidate["readId"] == "read1"
    assert candidate["origSeq"] == "AAAACCCCGGGG"
    assert candidate["seq"] == "AAAACCCCGGGG"
    assert candidate["strand"] == "fwd"
    assert candidate["kmer"] == "AAAA"
    assert candidate["window"] == 0
    assert candidate["k_start"] == 0
    assert candidate["k_end"] == 4
    assert candidate["read_len"] == 12

def test_pref_coords_window_zero():
    start, end = pref_coords(win=0, len_kmer=4)

    assert start == 0
    assert end == 4

def test_pref_coords_window_one():
    start, end = pref_coords(win=1, len_kmer=4)

    assert start == 1
    assert end == 5

def test_suf_coords_window_zero():
    seq = "AAAACCCCGGGG"

    start, end = suf_coords(win=0, len_kmer=4, seq=seq)

    assert start == 8
    assert end == 12

def test_suf_coords_window_one():
    seq = "AAAACCCCGGGG"

    start, end = suf_coords(win=1, len_kmer=4, seq=seq)

    assert start == 7
    assert end == 11

def test_kmer_tables_forward_prefix_and_suffix_window_zero():
    headers = {"AAAACCCCGGGG": "read1"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=4,n=1)

    pref_hit = prefix_table["AAAA"][0]
    suff_hit = suffix_table["GGGG"][0]

    assert pref_hit["readId"] == "read1"
    assert pref_hit["strand"] == "fwd"
    assert pref_hit["window"] == 0
    assert pref_hit["k_start"] == 0
    assert pref_hit["k_end"] == 4

    assert suff_hit["readId"] == "read1"
    assert suff_hit["strand"] == "fwd"
    assert suff_hit["window"] == 0
    assert suff_hit["k_start"] == 8
    assert suff_hit["k_end"] == 12

def test_kmer_tables_forward_prefix_and_suffix_window_one():
    headers = {"AAAACCCCGGGG": "read1"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=4,n=1)

    pref_hit = prefix_table["AAAC"][0]
    suff_hit = suffix_table["CGGG"][0]

    assert pref_hit["readId"] == "read1"
    assert pref_hit["strand"] == "fwd"
    assert pref_hit["window"] == 1
    assert pref_hit["k_start"] == 1
    assert pref_hit["k_end"] == 5

    assert suff_hit["readId"] == "read1"
    assert suff_hit["strand"] == "fwd"
    assert suff_hit["window"] == 1
    assert suff_hit["k_start"] == 7
    assert suff_hit["k_end"] == 11

def test_kmer_tables_include_reverse_complement_hits():
    headers = {"AAAACCCCGGGG": "read1"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=4,n=0)

    rc_pref_hit = prefix_table["CCCC"][0]
    rc_suff_hit = suffix_table["TTTT"][0]

    assert rc_pref_hit["readId"] == "read1"
    assert rc_pref_hit["strand"] == "rc"
    assert rc_pref_hit["seq"] == "CCCCGGGGTTTT"
    assert rc_pref_hit["window"] == 0

    assert rc_suff_hit["readId"] == "read1"
    assert rc_suff_hit["strand"] == "rc"
    assert rc_suff_hit["seq"] == "CCCCGGGGTTTT"
    assert rc_suff_hit["window"] == 0

def test_short_read_does_not_create_short_kmers():
    headers = {"ACG": "short_read"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=4,n=1)

    assert prefix_table == {}
    assert suffix_table == {}