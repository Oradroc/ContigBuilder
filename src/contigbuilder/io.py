def Fasta2Dict(filepath):
    '''Read in single line fasta file and return 
     seqs - dictionary to get seqs from ids of seqids:sequence
     headers - dictionary to get ids from seqs sequence:seqids'''
    seqs = {}
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    for i in range(len(lines)):
        if lines[i].startswith(">"):
            key = lines[i][1:]
            seqs[key] = lines[i+ 1]
    #if there is a duplicate sequence only one id is kept
    headers=dict(zip(seqs.values(),seqs.keys()))
    return seqs,headers

def writeFasta(contig, output_path, header="contig1"):
    '''Write a fasta file of the final contig'''
    with open(output_path, "w") as f:
        f.write(f">{header}\n")
        f.write(contig + "\n")

def writeAln(aln_df, output_path):
    '''Write final mapped alignment df file'''
    aln_df.to_csv(output_path, sep="\t", index=False)