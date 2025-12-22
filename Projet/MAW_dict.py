import numpy as np
from collections import defaultdict
from xopen import xopen
from readfa import readfq
import time

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

def canonicalKmers(read,k):
    #Return the set of kmers in read
    res=defaultdict(set)
    l = len(read)
    for i in range(l-k):
        res[read[i:i+k]].add(read[i-1])
    return res
    
def main(k_max,filename):
    #Output in "filename.tsv" all MAWs of size <= k_max
    output=open(filename+".tsv","w")
    i=0
    for name in reads:
        i+=1
        #print(name+"-->" + str(i)+'\n')
        read=reads[name]
        #kmers_dict_list = canonicalKmersList(read,k_max)
        kmers_dict_list=[canonicalKmers(read,k) for k in range(1,k_max)]
        for k in range(2,k_max):
            maw_list = set()
            kmers_dict = kmers_dict_list[k-1]
            kmers_dict_sup = kmers_dict_list[k-2]
            for x in kmers_dict:
                w_s_set = kmers_dict[x]
                w_p_set = kmers_dict_sup[x[:-1]]
                maw_list.update({l+x for l in w_p_set.difference(w_s_set)})
            maw_list = {canonicalSequence(x) for x in maw_list if reverseComplement(x) in maw_list}
            if maw_list != set():
                output.write(name+"\t"+str(k+1)+"\t"+(",".join(maw_list))+"\n")

def time_testing(k_max,filename):
    #Outputs in "filename.tsv" all MAWs of size <= k_max
    #Prints out progress with elated time
    start_time = time.time()
    output=open(filename+".tsv","w")
    i=0
    for name in reads:
        mid_time=time.time()
        i+=1
        print(name+"-->" + str(i)+ " with time " + str(mid_time-start_time)+'\n')
        read=reads[name]
        kmers_dict_list=[canonicalKmers(read,k) for k in range(1,k_max)]
        for k in range(2,k_max):
            maw_list = set()
            w_s_set = kmers_dict_list[k-1]
            w_p_set = kmers_dict_list[k-2]
            for w_s in w_s_set:
                w_s_letters = w_s_set[w_s]
                w_p_letters = w_p_set[w_s[:-1]]
                for l in w_p_letters.difference(w_s_letters):
                    maw_list.add(l+w_s)
            maw_list = {canonicalSequence(x) for x in maw_list if reverseComplement(x) in maw_list}
            if maw_list != set():
                output.write(name+"\t"+str(k+1)+"\t"+(",".join(maw_list))+"\n")
    stop_time=time.time()
    return(stop_time-start_time)

def alphabet_checking():
    alphabet = set()
    k=0
    for name in reads:
        k=k+1
        print(name+"-->" + str(k)+'\n')
        l=len(reads[name])
        for i in range(l):
            alphabet.add(reads[name][i])
    return alphabet
