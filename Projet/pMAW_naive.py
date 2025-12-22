#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq
import time
from itertools import groupby

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

def absentAmount(k_mers_dict_dict,x):
    #Returns how many dictionaries in k_mers_dict_dict x is absent from
    res= 0
    for name in k_mers_dict_dict:
        if x not in k_mers_dict_dict[name]:
            res+=1
    return res

def main(k_max,p,filenam):
    #Outputs to "filename.tsv" all the pMAWs of length <= k_max
    #Prints out progress as well as time elated
    start_time = time.time()
    candidates = set()
    res = set()
    k_mers_dict_dict = {}
    k_mers_total = set()
    i=0
    for name in reads:
        i+=1
        mid_time = time.time()
        print("Adding k-mers of read "+name+" --> "+str(i)+" with time "+str(mid_time-start_time)+"\n")
        for k in range(1,k_max):
            s=canonicalKmers(reads[name],k)
            k_mers_dict_dict.update({name:s})
            k_mers_total.update(s)
    size=len(reads)
    for x in k_mers_total:
        for a in Sigma:
            if absentAmount(k_mers_dict_dict,a+x)/size > p:
                candidates.add(a+x)
    for x in candidates:
        if (x[1:] not in candidates) and (x[:-1] not in candidates):
            res.add(a+x)
    stop_time = time.time()
    output_pmaws = groupby(res,len)
    output = open(filename+".tsv","w")
    for l,pmaw_group in output_pmaws:
        output.write(str(l)+"\t"+(",".join(pmaw_group))+"\n")
    return(stop_time-start_time)
            
