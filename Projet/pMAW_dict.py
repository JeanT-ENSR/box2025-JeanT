#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq
import time
import csv
from itertools import groupby

##### Test file #####
reads = {}

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

def absentAmount(x):
    #Returns how many dictionaries in k_mers_dict_dict x is absent from
    res= 0
    for name in reads:
        if x not in reads[name]:
            res+=1
    return res

visited = set()

def candidates(k_max,p,w):
    #Returns the set of potential pMAWs that are extensions of w
    global reads
    l=len(w)
    global visited
    s=len(reads)
    a=absentAmount(w)/s
    if l<= k_max and a>=p:
        return set(w)
    if l >= k_max:
        return set()
    c_set = set()
    for alpha in Sigma:
        x = alpha + w
        if not x in visited:
            b = absentAmount(x)/s
            if b>=p:
                c_set.add(x)
                visited.add(x)
            else:
                c_set.update(candidates(k_max,p,x))
        x = w + alpha
        if not x in visited:
            b = absentAmount(x)/s
            if b>=p:
                c_set.add(x)
                visited.add(x)
            else:
                c_set.update(candidates(k_max,p,x))
    return c_set

def get_maws(filename):
    maws_dict = defaultdict(set)
    with open(filename+".tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            maws_dict[row[0]].update(set(row[2].split(",")))
    return maws_dict

maws_dict = defaultdict(set)

def read_files(filename):
    global reads
    global maws_dict
    with xopen("all_ebi_plasmids.fa.xz") as fq:
        for name,seqfq,_ in readfq(fq):
            reads[name] = seqfq
    maws_dict = get_maws(filename)

def all_subwords(maw):
    #Yield all sub words of the input string
    for i in range(len(maw)):
        for j in range(i + 1, len(maw) + 1):
            if len(maw[i:j]) < len(maw):
                yield maw[i:j]

def main(k_max,p,filename):
    #Outputs to "filename.tsv" all the pMAWs of length <= k_max
    #Returns elated time
    start_time = time.time()
    c_set = set()
    res = set()
    for name in reads:
        for w in maws_dict[name]:
            c_set.update(candidates(k_max,p,w))
    for w in c_set:
        ispMAW = True
        for sub in all_subwords(w):
            if sub in c_seq:
                ispMAW = False
                break
        if ispMAW:
            res.add(w)
    stop_time = time.time()
    output_pmaws = groupby(res,len)
    output = open(filename+".tsv","w")
    for l,pmaw_group in output_pmaws:
        output.write(str(l)+"\t"+(",".join(pmaw_group))+"\n")
    return(stop_time-start_time)
            
