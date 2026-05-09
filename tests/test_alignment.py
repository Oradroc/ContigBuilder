import pytest
from contigbuilder.alignment import GetFinalCoords

def make_rec(read_id, seq, strand, qstart0, qend0, overlap_len=8, score=10):
    return {
        "hit": {
            "readId": read_id,
            "seq": seq,
            "strand": strand
        },
        "qstart0": qstart0,
        "qend0": qend0,
        "overlap_len": overlap_len,
        "score": score
    }

def test_right_extension_no_left_shift():
    Rthistory = [make_rec(read_id="read1",seq="CCCCGGGGTT",strand="fwd",qstart0=4,qend0=14)]
    Lefthistory = []

    df = GetFinalCoords(Rthistory, Lefthistory, start=0)
    row = df.iloc[0]

    assert row["sseqid"] == "read1"
    assert row["qseqid"] == "contig1"
    assert row["sstart"] == 1
    assert row["send"] == 10
    assert row["qstart"] == 5
    assert row["qend"] == 14
    assert row["strand"] == "fwd"
    assert row["overlap_len"] == 8
    assert row["score"] == 10

def test_left_extension_with_negative_start():
    Rthistory = []
    Lefthistory = [make_rec(read_id="read2",seq="GGGGAAAACC",strand="fwd",qstart0=-4,qend0=6)]

    df = GetFinalCoords(Rthistory, Lefthistory, start=-4)
    row = df.iloc[0]

    assert row["sseqid"] == "read2"
    assert row["sstart"] == 1
    assert row["send"] == 10
    assert row["qstart"] == 1
    assert row["qend"] == 10
    assert row["strand"] == "fwd"

def test_right_extension_after_left_shift():
    Rthistory = [make_rec(read_id="read3",seq="CCCCGGGGTT",strand="fwd",qstart0=4,qend0=14)]
    Lefthistory = []

    df = GetFinalCoords(Rthistory, Lefthistory, start=-4)
    row = df.iloc[0]

    assert row["qstart"] == 9
    assert row["qend"] == 18

def test_reverse_complement_read_coordinates():
    Rthistory = []
    Lefthistory = [make_rec(read_id="read4",seq="GGGGAAAACC",strand="rc",qstart0=-4,qend0=6)]

    df = GetFinalCoords(Rthistory, Lefthistory, start=-4)
    row = df.iloc[0]

    assert row["sseqid"] == "read4"
    assert row["sstart"] == 10
    assert row["send"] == 1
    assert row["qstart"] == 1
    assert row["qend"] == 10
    assert row["strand"] == "rc"

def test_combined_right_and_left_histories():
    Rthistory = [make_rec(read_id="right_read",seq="CCCCGGGGTT",strand="fwd",qstart0=4,qend0=14)]
    Lefthistory = [make_rec(read_id="left_read",seq="GGGGAAAACC",strand="fwd",qstart0=-4,qend0=6)]

    df = GetFinalCoords(Rthistory, Lefthistory, start=-4)

    assert len(df) == 2

    right_row = df[df["sseqid"] == "right_read"].iloc[0]
    left_row = df[df["sseqid"] == "left_read"].iloc[0]

    assert right_row["qstart"] == 9
    assert right_row["qend"] == 18

    assert left_row["qstart"] == 1
    assert left_row["qend"] == 10

def test_invalid_strand_raises_error():
    Rthistory = [make_rec(read_id="bad_read",seq="CCCCGGGGTT",strand="bad",qstart0=0,qend0=10)]
    Lefthistory = []

    with pytest.raises(ValueError):
        GetFinalCoords(Rthistory, Lefthistory, start=0)