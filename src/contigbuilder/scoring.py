from contigbuilder.sequenceUtils import GetMismatches

def bestseq(contig,hits,direction,len_kmer,mRate=0.1):
    '''Get best sequence from a list of candidates
    direction:
        right - means candidate extends the right side of the contig
        left - means candidate extends the left side of the contig
    Score:
    mRate-skip any candidates with a mismatch rate greater than this threshold start at 0.1
    len(overlap) * (1-mismatch_rate) + 0.5*additional sequence added
    Returns one best result dictionary, or None if no usable candidate exists.
    '''
    best = None

    for hit in hits:
        candidateSeq = hit["seq"]
        
        expectedOverlap = len_kmer+hit["window"]
        #build right
        if direction == "right":
            contig_overlap = contig[-expectedOverlap:]
            candidate_overlap = candidateSeq[:expectedOverlap]
            extension = candidateSeq[expectedOverlap:]
            new_contig = contig+extension
        #build left
        elif direction == "left":
                contig_overlap = contig[:expectedOverlap]
                candidate_overlap = candidateSeq[-expectedOverlap:]
                extension = candidateSeq[:-expectedOverlap]
                new_contig = extension+contig
        else:
            raise ValueError("direction must be right or left")
        mismatches= GetMismatches(contig_overlap, candidate_overlap)
        mismatch_rate = mismatches/expectedOverlap
        #skip candidates with poor mismatch
        if mismatch_rate>mRate:
            continue

        lengthAdded=len(extension)
        #if no length is added skip
        if lengthAdded==0:
            continue

        score= expectedOverlap*(1-mismatch_rate)+0.5*lengthAdded
        result ={
                "hit": hit,
                "seq": candidateSeq,
                "direction": direction,
                "overlap_len": expectedOverlap,
                "mismatches": mismatches,
                "mismatch_rate": mismatch_rate,
                "extension": extension,
                "lengthAdded": lengthAdded,
                "score": score,
                "new_contig": new_contig
                }
        if best is None or result['score']>best["score"]:
            best=result
    return best