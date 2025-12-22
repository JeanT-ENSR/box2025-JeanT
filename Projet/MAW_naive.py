#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq
import time

##### Test file #####
reads = {}
def read_files():
    global reads
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

##### MAW naive #####

def absentWord_aux(x,kmers):
    '''
    Assumption :
        k is egal to len(x)
        kmers is the dictionary of all kmer from read in canonical form.
    '''
    return not(canonicalSequence(x) in kmers)

def absentWord(read,x):
    return absentWord_aux(x,canonicalKmers(read,len(x)))

def minAbsentWord_aux(read,x,kmers,k1mers):
    assert(len(x) > 0)
    return absentWord_aux(x,kmers) and not(absentWord_aux(x[1:],k1mers)) and not(absentWord_aux(x[:-1],k1mers))

def minAbsentWord(read,x):
    l = len(x)
    assert(l > 0)
    kmers = canonicalKmers(read,l)
    k1mers = canonicalKmers(read,l-1)
    return minAbsentWord_aux(x,kmers,k1mers)

def allMinAbsentWord(read,kmax):
    '''
    2 steps :
        compute all kmers for all k less than kmax
        compute all MAW using the idea of the paper : "Scientific_Report.pdf"
    '''

    # Step 1
    kmers_array = []
    for k in range(kmax+1):
        kmers_array.append(canonicalKmers(read,k))
    
    # Step 2
    MAW = defaultdict(set)
    for k in range(1,kmax+1):
        for w in kmers_array[(k-1)]:
            for a in Sigma:
                s = w+a
                if(absentWord_aux(s,kmers_array[k]) and not(absentWord_aux(s[:-1], kmers_array[(k-1)]))):
                    MAW[k].add(s)
                s = a+w
                if(absentWord_aux(s,kmers_array[k]) and not(absentWord_aux(s[1:], kmers_array[(k-1)]))):
                    MAW[k].add(s)
    return MAW

##### test and output in a TSV #####
def main(kmax,filename)
    #Outputs to "filename.tsv" all the MAWs of length <= kmax
    #Prints progress as well as time elated
    start_time = time.time()
    aMAW_list = []
    sorted(reads.items(), key=lambda t: t[0])
    with open(filename+".tsv", "w", encoding="utf-8") as f:
        i=0
        for name in reads:
            i+=1
            mid_time = time.time()
            print(name+ " --> "+str(i)+" with time "+str(mid_time-start_time)+'\n')
            aMAW = allMinAbsentWord(reads[name], kmax)
            aMAW_list.append(aMAW)

            # Write in a TSV
            for k in aMAW:
                str_list = " ; ".join(aMAW[k])
                f.write(name + '\t' + str(k) + '\t' + str_list + "\n")
    stop_time = time.time()
    print("Job Done ! ")
    print(stop_time-start_time)

