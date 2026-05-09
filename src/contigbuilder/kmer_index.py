from tqdm import tqdm
from contigbuilder.sequenceUtils import prefix, suffix, RvsComp

def createCandidate(read,seq,kmer,win,headers,strand,k_start,k_end):
    '''
    create kmer candidate, take in orig rread r, seq, kmer, window, and header dict
    Return candidate with all items stored
    readId - id of the read
    origSeq- orignal read seq
    seq- oriented sequence orig read or rvs comp
    strand- "fwd" or "rc" (reverse comp)
    kmer-kmer sequence
    win-window offset of the kmer in the origSeq
    headers -dict mapping the orignal seq to the read id
    k_start,k_end- zero based, end exclusize coords of kmer in the oriented seq
    '''
    candidate={
    "readId": headers[read],
    "origSeq": read,
    "seq": seq,
    "strand": strand,
    "kmer": kmer,
    "window": win,
    "k_start": k_start,
    "k_end": k_end,
    "read_len": len(read)
    }
    return candidate
    
def pref_coords(win,len_kmer):
    '''return start and stop of prefix coords in kmer
    zero based end exclusive'''
    s=win
    e=win+len_kmer
    return s,e
    
def suf_coords(win,len_kmer,seq):
    '''return start and stop of suffix coords in kmer
    zero based end exclusive'''
    s=len(seq)-len_kmer-win 
    e=len(seq)-win
    return s,e
    
def get_NKmers_PrefSuftables(headers,len_kmer,n):
    '''
    return prefix table, suffix table
    Create both prefix and suffix lookup tables with all window kmers
    Get kmers for n steps.
    For example, n of 5 would get a 5 window kmer and all kmers up to the 5 nucleotide 
    '''
    prefix_table={}
    suffix_table={}

    reads=headers.keys()

    for read in tqdm(reads):
        #get both orientations
        orig=read
        oriented_reads = [("fwd", read),("rc", RvsComp(read))]
        for win in range(n+1):
            #based on fwd or rvs get pref and suffix and store in lookup table
            for orientation in oriented_reads:
                fwdrvs=orientation[0]#getr orientation fwd or rc
                r=orientation[1]#get the read
                #just in case if the read is shorter than the actual kmer skip it
                if len(r)<len_kmer+win:
                    continue
                #prefix
                pref=prefix(r,len_kmer,win)
                ps,pe = pref_coords(win,len_kmer)
                pref_candidate=createCandidate(orig,r,pref,win,headers,fwdrvs,ps,pe)
                if prefix_table.get(pref,None)==None:
                    prefix_table[pref]=[pref_candidate]
                else:
                    prefix_table[pref].append(pref_candidate)
                #suffix
                suf=suffix(r,len_kmer,win)
                ss,se = suf_coords(win,len_kmer,r)
                suf_candidate=createCandidate(orig,r,suf,win,headers,fwdrvs,ss,se)
                if suffix_table.get(suf,None)==None:
                    suffix_table[suf]=[suf_candidate]
                else:
                    suffix_table[suf].append(suf_candidate)
                
    return prefix_table,suffix_table