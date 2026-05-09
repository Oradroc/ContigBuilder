def prefix(s,len_kmer=5,n=0):
    '''
    Return the beginning of len elements of a string
    n - number of offset windoes 
    ex: ACGTCA for len 2 with n=0 would return AC n=1 CG
    '''
    return s[n:len_kmer+n]
def suffix(s,len_kmer=5,n=0):
    '''
    Return the end of len elements of a string
    n - number of offset windoes 
    ex: ACGTCA for len 2 with n=0 would return CA n=1 TC
    '''
    start=-(len_kmer+n)
    if n==0:
        return s[start:]
    else:
        return s[start:-n]
    
def RvsComp(s):
    '''
    Return the reverse complement of the string
    based on base pair rules
    '''
    pairs={"A":"T","T":"A","G":"C","C":"G"}
    rvs = s[::-1]
    newseq=""
    for n in rvs:
        newseq=newseq+pairs[n]
    return newseq

def GetMismatches(seq1, seq2):
    '''
    Return number of mismatches between two equal length strings
    Counts mismatches between strings.
    '''
    if len(seq1) != len(seq2):
        raise ValueError("seq1 and seq2 must be the same length")
    mismatches=0
    for a, b in zip(seq1, seq2):
        if a != b:
            mismatches += 1
    return mismatches