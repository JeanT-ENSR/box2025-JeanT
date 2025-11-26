#! /usr/bin/env python3

from xopen import xopen
from readfa import readfq
import numpy as np
# from collections import defaultdict

##### Utility function #####

def change(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    return c

def reverseComplement(read):
    s = list(read)
    l = len(s)
    for i in range(l):
        s[i] = change(s[i])
    return ''.join(list(reversed(s))) 

def open(fastaFilePath):
    sequences = []
    with xopen(fastaFilePath) as fa:
        for _,seq,_ in readfq(fa):
            sequences.append(seqfa)
    return sequences


##### MAW function #####

def isXInSeq(x,seq):
    '''
    Return True iff x or its reverse complement is in seq
    '''
    isIn = False
    k = len(x)
    l = len(seq)
    reverseX = reverseComplement(x)
    for i in range(l-k):
        subSeq = seq[i:i+k]
        isIn = isIn or (subSeq == x) or (subSeq == reverseX)
    return isIn

# Test passed
assert(isXInSeq("AC","TTACTT"))
assert(isXInSeq("AC","GTTG"))
assert(not isXInSeq("AC","TTTTTT"))
assert(not isXInSeq("AC","GATCG"))

def subSeq(sequence,kmin,kmax):
    l = len(sequence)
    subSequence = {}
    for k in range(kmin,kmax):
        for i in range(l-k):
            subSequence[sequence[i:i+k]] = True
    return subSequence

def subSeqSequences(sequences,kmin,kmax):
    subSeqSequences = defaultdict(list) # Need to store the name of the sequences
    for sequence in sequences:
        subSequence = subSeq(sequence,kmin,kmax)
        for sub in subSequence:
            subSeqSequences[sub].append(sub)

