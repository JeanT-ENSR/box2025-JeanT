#! /usr/bin/env python3

from xopen import xopen
from readfa import readfq
import numpy as np
import sys
from collections import defaultdict

new_limit = 4 * 40000 * 150
sys.setrecursionlimit(new_limit)  

reads = []
refs = []

with xopen("reads.fq.gz") as fq:
    for _,seqfq,_ in readfq(fq):
        reads.append(seqfq)


with xopen("ref.fa.gz") as fa:
    for _,seqfa,_ in readfq(fa):
        refs.append(seqfa)

##### Read Aligner NAIF #####

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
    
    operations = []
    for path in paths:
        (j0,i0) = path[0]
        operation = []
        l = len(path)
        for i in range(1,l):
            (j,i) = path[i]
            (di,dj) = (i0-i,j0-j)
            if (di,dj) == (1,1):
                if memo[i0][j0] == memo[i][j] :
                    operation.append('=')
                else:
                    operation.append('X')
            elif (di,dj) == (1,0):
                operation.append('I')
            else:
                operation.append('D')
            (j0,i0) = (j,i)
        operations.append(operation)

    return operations
            
# p = "TACGTCAGT"
# t = "AACCCTATGTCATGCCTTGGA"
# a = approximateMatching(p,t)
# print(a)
# print(buildEditOperation(a))

##### Read Aligner Seed-and-Extend #####

k = 29
kMers = defaultdict(list)
nb_refs = len(refs)
for j in range(nb_refs):
    l = len(refs[j])
    for i in range(l):
        if i+k < l:
            kMers[refs[j][i:i+k]].append((j,i))

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

# read_test = "ATTTACCGTA"
# print(read_test)
# print(reverseComplement(read_test))

def seed(read):
    seed = []
    l = len(read)
    for i in range(l):
        if i+k < l:
            if not(not(kMers[read[i:i+k]])):
                seed.append((kMers[read[i:i+k]],i,True)) # True = sens d'origine du read
    reverseRead = reverseComplement(read)
    for i in range(l):
        if i+k < l:
            if not(not(kMers[reverseRead[i:i+k]])):
                seed.append((kMers[reverseRead[i:i+k]],i,False))  # False = reverse du read
    return seed

def extend(seed,read):
    '''
    Si seed est vide alors on passe et on donne d = len(read) et full deletion
    Si seed n'est pas vide alors
    Je le fais en 2 étapes :
        - Je cherche le plus grand perfect match
        - Je fais un approximativeMatch entre le read et une zone autour du perfect match dans la ref
    '''
    if len(seed) == 0 :
        l = len(read)
        return ['X' for _ in range(l)]
    # On cherche le plus grand perfect match grace au seed
    # Idée : Si deux kmer coté à coté sont dans le seed alors on a un grand kmer de taille k+1
    newSeed = []
    for position in seed:
        (refIndexList, indexRead, direction) = position
        for (numRef, refIndex) in refIndexList:
            isIn = False
            newSeed2 = []
            for newPosition in newSeed:
                (numRef2, refStart, refIndex2, indexRead2, direction2) = newPosition
                if (direction == direction2) and (numRef == numRef2) and (refIndex2 + 1 == refIndex) and (indexRead2 + 1 == indexRead):
                    newSeed2.append((numRef, refStart, refIndex, indexRead, direction))
                    isIn = True
                else:
                    newSeed2.append(newPosition)
            newSeed = newSeed2
            if not(isIn):
                newSeed.append((numRef, refIndex, refIndex, indexRead, direction))
    
    # Etape 2
    for (numRef, refIndexStart, refIndex, readIndex, direction) in newSeed:
        refIndexEnd = refIndex + k
        readIndexEnd = readIndex + k
        bigK = (refIndexEnd - refIndexStart)
        readIndexStart = readIndexEnd - bigK
        l = len(read)
        lRef = len(refs[numRef])
        refStart = max(refIndexStart - (2 * readIndexStart), 0)
        refEnd = min(refIndexEnd + (2 * (l - readIndexEnd)),lRef)

        if direction :
            memo = approximateMatching(read,refs[numRef][refStart:refEnd])
            return buildEditOperation(memo)
        else:
            memo = approximateMatching(reverseComplement(read),refs[numRef][refStart:refEnd])
            return buildEditOperation(memo)
        

for read in reads:
    t = extend(seed(read),read)
    print(t)
