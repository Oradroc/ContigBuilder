from contigbuilder.sequenceUtils import prefix, suffix
from contigbuilder.scoring import bestseq
from contigbuilder.kmer_index import get_NKmers_PrefSuftables

def buildRight(prefix_table,Query,len_kmer,mRate=0.1,maxsteps=10000):
    '''
    return contig built right with addition history, start and stop of contig
    Using prefix table build right on the suffix of the query'''
    contig = Query
    history = []
    used_reads = set()
    #get initial coordinates for contig
    StartContig=0
    EndContig=len(Query)
    
    for step in range(maxsteps):
        QSuf = suffix(contig, len_kmer, 0)
        #get kmer if no kmer is found return empty list
        hits = prefix_table.get(QSuf, [])
        #Remove any reads that have already been used
        unusedHits = []
        for hit in hits:
            if hit["readId"] in used_reads:
                continue
            unusedHits.append(hit)
        #Move on if all the reads for the kmer were used and return the final contig
        if len(unusedHits) == 0:
            break
        #get the best sequence for list of hits
        best = bestseq(contig=contig,hits=unusedHits,direction="right",len_kmer=len_kmer,mRate=mRate)
        #if no best hit is found because the remaining hits have bad overlap quality, then break
        if best is None:
            break
        
        read_len=len(best["hit"]["seq"])
        #get start of the initial overlap
        qstart0=EndContig-best["overlap_len"]
        #get end
        qend0=qstart0+read_len
        best["qstart0"]=qstart0
        best["qend0"]=qend0
        #update end of contig start stays the same
        EndContig=qend0

        contig = best["new_contig"]

        used_reads.add(best["hit"]["readId"])
        history.append(best)

    return contig, history, StartContig,EndContig,used_reads

def buildLeft(suffix_table, Query, len_kmer,used_reads=None,mRate=0.1, maxsteps=10000):
    '''
    Take in the output of the build right contig and build left 
    return the final contig, history, and start and stop
    Using suffix table build left on the prefix of the query'''
    contig = Query
    history = []
    if used_reads is None:
        used_reads = set()
    #get initial coordinates for contig
    StartContig=0
    EndContig=len(Query) #end of right build

    for step in range(maxsteps):
        QPref = prefix(contig, len_kmer, 0)
        #get kmer if no kmer is found return empty list
        hits = suffix_table.get(QPref, [])
        #Remove any reads that have already been used
        unusedHits = []
        for hit in hits:
            if hit["readId"] in used_reads:
                continue
            unusedHits.append(hit)
        #Move on if all the reads for the kmer were used and return the final contig
        if len(unusedHits) == 0:
            break
        #get the best sequence for list of hits
        best = bestseq(contig=contig,hits=unusedHits,direction="left",len_kmer=len_kmer,mRate=mRate)
        #if no best hit is found because the remaining hits have bad overlap quality, then break
        if best is None:
            break

        read_len=len(best["hit"]["seq"])
        #add the overlap len to get the neg end
        qend0=StartContig+best["overlap_len"]
        #subtract length from the end to get moving edn
        qstart0=qend0-read_len

        best["qstart0"]=qstart0
        best["qend0"]=qend0

        contig=best["new_contig"]
        StartContig=qstart0

        used_reads.add(best["hit"]["readId"])
        history.append(best)
    return contig, history, StartContig,EndContig

def BuildContig(headers, Query, k_right=7, n_right=6, k_left=7, n_left=5, mRate=0.1, maxsteps=10000):
    '''
    Build final contig by extending right first, then left.
    Allows different kmer/window settings for right and left extension.
    '''
    #if both are the same only need to call function once
    if (k_right==k_left) and (n_right==n_left):
        prefix_table, suffix_table = get_NKmers_PrefSuftables(headers, k_right, n_right)
    #if different, then need to call separate
    else:
        prefix_table, _ = get_NKmers_PrefSuftables(headers, k_right, n_right)
        _, suffix_table = get_NKmers_PrefSuftables(headers, k_left, n_left)

    Rtcontig, Rthistory, Rstart, Rend, used_reads = buildRight(prefix_table=prefix_table,Query=Query,len_kmer=k_right,mRate=mRate,maxsteps=maxsteps)

    FinalContig, Lefthistory, Lstart, Lend = buildLeft(suffix_table=suffix_table,Query=Rtcontig,len_kmer=k_left,used_reads=used_reads,mRate=mRate,maxsteps=maxsteps)

    return FinalContig, Rthistory, Lefthistory, Rstart, Rend, Lstart, Lend