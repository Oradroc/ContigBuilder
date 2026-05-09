import pandas as pd
def GetFinalCoords(Rthistory,Lefthistory,start):
    '''Convert internally tracked zero based, end exclusive contig coordinates
    into one based inclusive coordinates for the ALLELES.aln file.'''
    records=Rthistory+Lefthistory
    shift=-start #get the shift from the neg start of the left buildLeft
    rows=[]
    for rec in records:
        hit=rec["hit"]
        qstart=rec["qstart0"]+shift+1
        qend=rec["qend0"]+shift
        read_len=len(hit["seq"])

        #get seq coords for if fwd or revs comp
        if hit["strand"]=="fwd":
            sstart=1
            send=read_len
        #make star opposite to show original strand coords
        elif hit["strand"]=="rc":
            sstart=read_len
            send=1
        else:
            raise ValueError("strand must be fwd or rc")
        
        row = {
            "sseqid": hit["readId"],
            "qseqid": "contig1",
            "sstart": sstart,
            "send": send,
            "qstart": qstart,
            "qend": qend,
            "strand": hit["strand"],
            "overlap_len": rec["overlap_len"],
            "score": rec["score"]
        }
        rows.append(row)
    return pd.DataFrame(rows)