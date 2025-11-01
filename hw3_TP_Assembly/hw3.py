#! /usr/bin/env python3

from readfa import readfq
from xopen import xopen

def create_dgb(fasta_file, k, t):
    print("Je n'ai pas enlever le reverse complement")

    file = open(fasta_file, 'r')
    reads = list(readfq(file))
    seqs = []
    for i in range(len(reads)):
        seqs.append(reads[i][1])

    # creation du dictionnnaire pour chaque k-mer
    k_mer = {}
    for seq in seqs:
        length_seq = len(seq)
        for i in range(length_seq-k+1):
            buffer = seq[i:i+k]
            if (buffer in k_mer):
                k_mer[buffer] += 1
            else:
                k_mer[buffer] = 1

    # nettoyage du dictionnaire avec t
    k_mer_list = []
    for s in k_mer.keys():
        k_mer_list.append(s)
    for s in k_mer_list:
        if k_mer[s] < t:
            del k_mer[s]

    return k_mer

### ecoli_genome_150k ###
print("ecoli_genome_150k")

k_mers_ecoli_genome = create_dgb("ecoli_genome_150k.fa", 31, 1)
print("nb k-mer stored in the graph : ", end='')
print(len(k_mers_ecoli_genome))

nb_processed = 0
for k in k_mers_ecoli_genome:
    nb_processed += k_mers_ecoli_genome[k]
print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
print(nb_processed)


k_min = 16
k_max = 17
t = 3

print("\nt = ", end='')
print(t)
### ecoli_sample_perfect_reads_forward ###
print("\necoli_sample_perfect_reads_forward")

for k in range(k_min,k_max):
    print("k = ", end='')
    print(k)

    k_mers_ecoli_sample_perfect_reads_forward = create_dgb("reads/ecoli_sample_perfect_reads_forward.fasta", k, t)
    print("nb k-mer stored in the graph : ", end='')
    print(len(k_mers_ecoli_sample_perfect_reads_forward))

    nb_processed = 0
    for key in k_mers_ecoli_sample_perfect_reads_forward:
        nb_processed += k_mers_ecoli_sample_perfect_reads_forward[key]
    print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
    print(nb_processed)


### ecoli_sample_perfect_reads ###
print("\necoli_sample_perfect_reads")

for k in range(k_min,k_max):
    print("k = ", end='')
    print(k)

    k_mers_ecoli_sample_perfect_reads = create_dgb("reads/ecoli_sample_perfect_reads.fasta", k, t)
    print("nb k-mer stored in the graph : ", end='')
    print(len(k_mers_ecoli_sample_perfect_reads))

    nb_processed = 0
    for key in k_mers_ecoli_sample_perfect_reads:
        nb_processed += k_mers_ecoli_sample_perfect_reads[key]
    print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
    print(nb_processed)


### ecoli_sample_reads_01 ###
print("\necoli_sample_reads_01")

for k in range(k_min,k_max):
    print("k = ", end='')
    print(k)

    k_mers_ecoli_sample_reads_01 = create_dgb("reads/ecoli_sample_reads_01.fasta", k, t)
    print("nb k-mer stored in the graph : ", end='')
    print(len(k_mers_ecoli_sample_reads_01))

    nb_processed = 0
    for key in k_mers_ecoli_sample_reads_01:
        nb_processed += k_mers_ecoli_sample_reads_01[key]
    print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
    print(nb_processed)


### ecoli_sample_reads ###
print("\necoli_sample_reads")

for k in range(k_min,k_max):
    print("k = ", end='')
    print(k)

    k_mers_ecoli_sample_reads = create_dgb("reads/ecoli_sample_reads.fasta", k, t)
    print("nb k-mer stored in the graph : ", end='')
    print(len(k_mers_ecoli_sample_reads))

    nb_processed = 0
    for key in k_mers_ecoli_sample_reads:
        nb_processed += k_mers_ecoli_sample_reads[key]
    print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
    print(nb_processed)