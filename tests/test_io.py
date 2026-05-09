import pandas as pd

from contigbuilder.io import Fasta2Dict, writeFasta, writeAln

def test_Fasta2Dict_reads_single_sequence(tmp_path):
    fasta_path = tmp_path / "test.fasta"

    fasta_path.write_text(">seq1\nACGTCA\n")

    seqs, headers = Fasta2Dict(fasta_path)

    assert seqs == {"seq1": "ACGTCA"}
    assert headers == {"ACGTCA": "seq1"}

def test_Fasta2Dict_reads_multiple_sequences(tmp_path):
    fasta_path = tmp_path / "test.fasta"

    fasta_path.write_text(">read1\nACGT\n>read2\nTTTT\n")

    seqs, headers = Fasta2Dict(fasta_path)

    assert seqs["read1"] == "ACGT"
    assert seqs["read2"] == "TTTT"

    assert headers["ACGT"] == "read1"
    assert headers["TTTT"] == "read2"

def test_writeFasta_writes_expected_file(tmp_path):
    out_path = tmp_path / "ALLELES.fasta"

    writeFasta(contig="ACGTCATGCA",output_path=out_path,header="contig1")

    text = out_path.read_text()

    assert text == ">contig1\nACGTCATGCA\n"

def test_writeAln_writes_tab_separated_file(tmp_path):
    out_path = tmp_path / "ALLELES.aln"

    aln_df = pd.DataFrame([
        {
            "sseqid": "read1",
            "qseqid": "contig1",
            "sstart": 1,
            "send": 10,
            "qstart": 5,
            "qend": 14
        }
    ])

    writeAln(aln_df, out_path)

    text = out_path.read_text()

    assert "sseqid\tqseqid\tsstart\tsend\tqstart\tqend" in text
    assert "read1\tcontig1\t1\t10\t5\t14" in text