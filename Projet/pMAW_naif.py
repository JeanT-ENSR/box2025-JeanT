#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq

##### Test file #####
reads = {}
with xopen("all_ebi_plasmids.fa.xz") as fq:
    for name,seqfq,_ in readfq(fq):
        reads[name] = seqfq

##### Utils #####

Sigma = {'A', 'C', 'G', 'T'}

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

def canonicalSequence(read):
    return min(read,(reverseComplement(read)))

def canonicalKmers(read,k):
    kmers = set()
    l = len(read)
    for i in range(l):
        if i+k-1 < l:
            kmers.add(canonicalSequence(read[i:i+k]))
    return kmers

def absentWord_aux(x,kmers):
    '''
    Assumption :
        k is egal to len(x)
        kmers is the dictionary of all kmer from read in canonical form.
    '''
    return not(canonicalSequence(x) in kmers)

##### PMAW naive #####
def allPMinAbsentWord(reads,p,kmax):
    assert(0 < p)
    assert(p < 1)
    nb = int(p*len(reads))

    kmersReads_array = defaultdict(list)
    for name in reads:
        print(name)
        for k in range(kmax+1):
            kmersReads_array[name].append(canonicalKmers(reads[name],k))
    print("kmers done")

    pMAWReads = defaultdict(defaultdict(int))
    for name in reads:
        for k in range(1,kmax+1):
            for w in kmersReads_array[name][(k-1)]:
                for a in Sigma:
                    isin = False
                    s = w+a
                    if(absentWord_aux(s,kmersReads_array[name][k]) and not(absentWord_aux(s[:-1], kmersReads_array[name][(k-1)]))):
                        isin = True
                    s = a+w
                    if(absentWord_aux(s,kmersReads_array[name][k]) and not(absentWord_aux(s[1:], kmersReads_array[name][(k-1)]))):
                        isin = True
                    if isin:
                        pMAWReads[k][s] = pMAWReads[k][s]+1

    filter_pMAWReads = defaultdict(int, {k: v for k, v in pMAWReads.items() if v > nb})
    return filter_pMAWReads

##### test and output in a TSV #####
kmax = 3
p = 0.5
sorted(reads.items(), key=lambda t: t[0])
with open("pMAW.tsv", "w", encoding="utf-8") as f:
    aPMAW = allPMinAbsentWord(reads, p, kmax)

    # Write in a TSV
    for k in aPMAW:
        str_list = " ; ".join(aPMAW[k])
        f.write(str(k) + '\t' + str_list + "\n")

print("Job Done ! ")

