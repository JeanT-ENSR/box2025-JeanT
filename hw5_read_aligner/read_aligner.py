#! /usr/bin/env python3

from xopen import xopen
from readfa import readfq
import numpy as np
import sys  

new_limit = 4 * 40000 * 150
sys.setrecursionlimit(new_limit)  

##### Read Aligner NAIF #####

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

def findPath(iStart,memo):
    path = []
    i = iStart
    j = len(memo) - 1
    path.append((i,j))
    while(j > 0):
        if i > 0:
            vMin = min(memo[j-1][i-1], memo[j][i-1], memo[j-1][i])
            if memo[j-1][i-1] == vMin:
                j = j-1
                i = i-1
                path.append((i,j))
            elif memo[j][i-1] == vMin:
                i = i-1
                path.append((i,j))
            else:
                j = j-1
                path.append((i,j))
        else:
            j = j-1
            path.append((i,j))
    return path


def buildEditOperation(memo):
    l1 = len(memo)
    l2 = len(memo[l1-1])
    dMin = min(memo[l1-1])
    paths = []
    for i in range(l2):
        if memo[l1-1][i] == dMin:
            paths.append(findPath(i,memo))
    return paths
            
p = "TACGTCAGT"
t = "AACCCTATGTCATGCCTTGGA"
a = approximateMatching(p,t)
print(a)
print(buildEditOperation(a))

##### Read Aligner Seed-and-Extend #####
