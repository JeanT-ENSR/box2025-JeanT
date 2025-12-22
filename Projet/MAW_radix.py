import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq
import time
import os.path
from radix_trees import *


##### Test file #####
name_test=""
reads = {}
def read_files():
    global name_test
    with xopen("all_ebi_plasmids.fa.xz") as fq:
        for name,seqfq,_ in readfq(fq):
            name_test=name
            reads[name] = seqfq

##### Utils #####

Sigma = {'A', 'C', 'G', 'T'}
alphabet = ['A','T','C','G']

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

def construct_tree(read,k_max):
    #Construct a radix tree containing every k_max-mer
    t=RadixTree()
    i=0
    l=len(read)
    while i+k_max<l: #Get as many k_max-mers as possible
        seq = read[i:i+k_max]
        t.insert_node(seq,read[i-1])
        i+=1
    while i<l: #We get the rest as k-mers where k<k_max
        seq = read[i:]
        t.insert_node(seq,read[i-1])
        i+=1
    return t

def getMAWList(t,k_max,l):
    #t is the tree to explore, k_max is the max length of an exploration, l is the list of nodes explored up until now
    #We compute a depth-first-search of the tree like in algorithm 5 in the report
    res = [set() for i in range(k_max)]
    l.append(t)
    current_node = t
    if t._next != {}:
        for k in t._next:
            c = t._next[k]
            temp = getMAWList(c,k_max,l)
            for i in range(k_max):
                res[i].update(temp[i])
    else:
        size=l[-1]._key_size
        k=1
        temp = l[0]
        for i in range(size):
            if i> temp._key_size:
                w_p_set = temp._data._data
                temp = l[k]
                k+=1
                w_s_set = temp._data._data
                letters_candidates = w_p_set.difference(w_s_set)
                for s in letters_candidates:
                    res[i].add(s+temp._key[:i])
    return res
            
        

def main(k_max,filename):
    #Outputs to "filename.tsv" all MAWs of length <= k_max
    output = open(filename+".tsv","w")
    i=0
    start_time = time.time()
    for name in reads:
        i+=1
        mid_time = time.time()
        print(name+"-->" + str(i)+" with time "+str(mid_time-start_time)+'\n')
        read=reads[name]
        t=construct_tree(read,k_max)
        res = getMAWList(t._tree,k_max,[])
        for k in range(k_max):
            res = {canonicalSequence(x) for x in res if reverseComplement(x) in res}
            if res != set():
                output.write(name+"\t"+str(k+1)+"\t"+(",".join(res))+"\n")
    end_time = time.time()
    return(end_time - start_time)
