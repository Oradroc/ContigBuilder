from contigbuilder.io import Fasta2Dict, writeFasta, writeAln
from contigbuilder.assembler import BuildContig
from contigbuilder.alignment import GetFinalCoords

def test_end_to_end_small_dataset(tmp_path):
    query_path = tmp_path / "QUERY.fasta"
    reads_path = tmp_path / "READS.fasta"
    out_fasta = tmp_path / "ALLELES.fasta"
    out_aln = tmp_path / "ALLELES.aln"

    query_path.write_text(">query1\nCCCCGGGG\n")

    reads_path.write_text(">left_read\nAAAACCCCGGGG\n"">right_read\nCCCCGGGGTTTT\n")

    QueryDict, _ = Fasta2Dict(query_path)
    Query = list(QueryDict.values())[0]

    Reads, headers = Fasta2Dict(reads_path)

    FinalContig, Rthistory, Lefthistory, Rstart, Rend, Lstart, Lend = BuildContig(headers=headers,Query=Query,k_right=8,n_right=0,k_left=8,n_left=0,mRate=0.1,maxsteps=10)

    aln_df = GetFinalCoords(Rthistory=Rthistory,Lefthistory=Lefthistory,start=Lstart)

    writeFasta(FinalContig, out_fasta, header="contig1")
    writeAln(aln_df, out_aln)

    assert FinalContig == "AAAACCCCGGGGTTTT"

    fasta_text = out_fasta.read_text()
    assert fasta_text == ">contig1\nAAAACCCCGGGGTTTT\n"

    aln_text = out_aln.read_text()
    assert "sseqid\tqseqid\tsstart\tsend\tqstart\tqend" in aln_text
    assert "left_read" in aln_text
    assert "right_read" in aln_text

    assert len(aln_df) == 2