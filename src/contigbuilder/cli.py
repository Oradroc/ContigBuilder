import argparse
from contigbuilder.io import Fasta2Dict, writeFasta, writeAln
from contigbuilder.assembler import BuildContig
from contigbuilder.alignment import GetFinalCoords

def parse_args():
    parser = argparse.ArgumentParser(description="Query seeded greedy contig builder")
    parser.add_argument("--query",required=True,help="Path to QUERY.fasta")
    parser.add_argument("--reads",required=True,help="Path to READS.fasta")
    parser.add_argument("--out_fasta",default="ALLELES.fasta",help="Output FASTA path")
    parser.add_argument("--out_aln",default="ALLELES.aln",help="Output alignment table path")
    parser.add_argument("--k_right",type=int,default=7,help="Kmer length for right extension")
    parser.add_argument("--n_right",type=int,default=6,help="Maximum window offset for right extension")
    parser.add_argument("--k_left",type=int,default=7,help="Kmer length for left extension")
    parser.add_argument("--n_left",type=int,default=5,help="Maximum window offset for left extension")
    parser.add_argument("--mRate",type=float,default=0.1,help="Maximum allowed mismatch rate")
    parser.add_argument("--maxsteps",type=int,default=10000,help="Maximum extension steps per direction")
    return parser.parse_args()

def main():
    args = parse_args()
    QueryDict, QueryHeader = Fasta2Dict(args.query)
    Query = list(QueryDict.values())[0]
    Reads, headers = Fasta2Dict(args.reads)
    FinalContig, Rthistory, Lefthistory, Rstart, Rend, Lstart, Lend = BuildContig(headers=headers,Query=Query, \
        k_right=args.k_right,n_right=args.n_right,k_left=args.k_left,n_left=args.n_left,mRate=args.mRate,maxsteps=args.maxsteps)
    aln_df = GetFinalCoords(Rthistory=Rthistory,Lefthistory=Lefthistory,start=Lstart)
    writeFasta(contig=FinalContig,output_path=args.out_fasta,header="contig1")
    writeAln(aln_df=aln_df,output_path=args.out_aln)
    print("Finished contig assembly")
    print("Final contig length:", len(FinalContig))
    print("Right extensions:", len(Rthistory))
    print("Left extensions:", len(Lefthistory))
    print("Wrote FASTA:", args.out_fasta)
    print("Wrote alignment:", args.out_aln)

if __name__ == "__main__":
    main()