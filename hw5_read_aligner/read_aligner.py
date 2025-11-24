#! /usr/bin/env python3

from xopen import xopen
from readfa import readfq
import numpy as np
import sys  

new_limit = 4 * 40000 * 150
sys.setrecursionlimit(new_limit)  

##### Read Aligner

reads = []
ref = []

with xopen("reads.fq.gz") as fq:
    for _,seqfq,_ in readfq(fq):
        reads.append(seqfq)


with xopen("ref.fa.gz") as fa:
    for _,seqfa,_ in readfq(fa):
        ref.append(seqfa)

max_int = max([len(read) for read in reads]) + 1

def approximateMatching(read,ref):
    p = len(read) + 1
    t = len(ref) + 1
    memo = np.full((p,t), max_int, dtype=np.int64)
    for i in range(t):
        memo[0][i] = 0
    for i in range(p):
        memo[i][0] = i
    edDistRecursiveMemo(read,ref,memo)
    return memo

def edDistRecursiveMemo(a, b, memo):
    if memo[len(a)][len(b)] < max_int :
        return memo[(len(a), len(b))]
    delt = 1 if a[-1] != b[-1] else 0
    ans = min(edDistRecursiveMemo(a[:-1], b[:-1], memo) + delt, edDistRecursiveMemo(a[:-1], b, memo) + 1, edDistRecursiveMemo(a, b[:-1], memo) + 1)
    memo[(len(a), len(b))] = ans
    return ans

def buildEditOperation(memo):
    return "TODO"

print(len(reads))
print(len(ref))

print(reads[0][:5])
print(reads[1][:12])

print(approximateMatching(reads[0][:5], reads[1][:12]))

p = "TACGTCAGT"
t = "AACCCTATGTCATGCCTTGGA"
print(approximateMatching(p,t))

match1 = approximateMatching(reads[0],ref[0])
match2 = approximateMatching(reads[0],ref[1])
print(min(match1[150]))
print(min(match2[150]))