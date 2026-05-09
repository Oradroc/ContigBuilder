from contigbuilder.assembler import buildRight, buildLeft, BuildContig
from contigbuilder.kmer_index import get_NKmers_PrefSuftables

def test_buildRight_single_extension():
    Query = "AAAACCCCGGGG"

    headers = {"CCCCGGGGTTTT": "read1"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=8,n=0)

    contig, history, start, end, used_reads = buildRight(prefix_table=prefix_table,Query=Query,len_kmer=8,mRate=0.1,maxsteps=10)

    assert contig == "AAAACCCCGGGGTTTT"
    assert len(history) == 1
    assert start == 0
    assert end == 16
    assert "read1" in used_reads

    best = history[0]
    assert best["hit"]["readId"] == "read1"
    assert best["qstart0"] == 4
    assert best["qend0"] == 16
    assert best["overlap_len"] == 8

def test_buildLeft_single_extension():
    Query = "CCCCGGGGTTTT"

    headers = {"AAAACCCCGGGG": "read1"}

    prefix_table, suffix_table = get_NKmers_PrefSuftables(headers=headers,len_kmer=8,n=0)

    contig, history, start, end = buildLeft(suffix_table=suffix_table,Query=Query,len_kmer=8,used_reads=set(),mRate=0.1,maxsteps=10)

    assert contig == "AAAACCCCGGGGTTTT"
    assert len(history) == 1
    assert start == -4
    assert end == 12

    best = history[0]
    assert best["hit"]["readId"] == "read1"
    assert best["qstart0"] == -4
    assert best["qend0"] == 8
    assert best["overlap_len"] == 8

def test_buildRight_no_hits_returns_original_query():
    Query = "AAAACCCCGGGG"
    prefix_table = {}

    contig, history, start, end, used_reads = buildRight(prefix_table=prefix_table,Query=Query,len_kmer=8,mRate=0.1,maxsteps=10)

    assert contig == Query
    assert history == []
    assert start == 0
    assert end == len(Query)
    assert used_reads == set()

def test_buildLeft_no_hits_returns_original_query():
    Query = "AAAACCCCGGGG"
    suffix_table = {}

    contig, history, start, end = buildLeft(suffix_table=suffix_table,Query=Query,len_kmer=8,used_reads=set(),mRate=0.1,maxsteps=10)

    assert contig == Query
    assert history == []
    assert start == 0
    assert end == len(Query)

def test_BuildContig_right_then_left():
    Query = "CCCCGGGG"

    headers = {"CCCCGGGGTTTT": "right_read","AAAACCCCGGGG": "left_read"}

    FinalContig, Rthistory, Lefthistory, Rstart, Rend, Lstart, Lend = BuildContig(headers=headers,Query=Query,k_right=8,n_right=0,k_left=8,n_left=0,mRate=0.1,maxsteps=10)

    assert FinalContig == "AAAACCCCGGGGTTTT"
    assert len(Rthistory) == 1
    assert len(Lefthistory) == 1

    assert Rstart == 0
    assert Rend == 12

    assert Lstart == -4
    assert Lend == 12