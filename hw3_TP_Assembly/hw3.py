#! /usr/bin/env python3

from readfa import readfq
# from xopen import xopen
import matplotlib.pyplot as plt
import numpy as np

print("Ma structure pour enregistrer un dgb est dictionnaire qui stocke le nombre de fois que chaque k-mer est vu dans le fichier\n")

def create_dgb(fasta_file, k, t):
    # print("Je n'ai pas enlever le reverse complement")

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
print("Q3 : ecoli_genome_150k")

k_mers_ecoli_genome = create_dgb("ecoli_genome_150k.fa", 31, 1)
print("nb k-mer stored in the graph : ", end='')
print(len(k_mers_ecoli_genome))

nb_processed = 0
for k in k_mers_ecoli_genome:
    nb_processed += k_mers_ecoli_genome[k]
print("nb k-mer processed is more than (because I clear some k_mers): ", end='')
print(nb_processed)


#### Q4 ####
print("\n\nQ4 :")
k_min = 30
k_max = 31
t = 0
t_max = 20

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
    print("nb k-mer processed is more than : ", end='')
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
    print("nb k-mer processed is more than : ", end='')
    print(nb_processed)

### ecoli_sample_reads_01 ###
print("\n\necoli_sample_reads_01 : ")
stored1 = np.zeros(t_max)
processed1 = np.zeros(t_max)
for t in range(t_max):
    print("\nt = ", end='')
    print(t)

    for k in range(k_min,k_max):
        print("k = ", end='')
        print(k)

        k_mers_ecoli_sample_reads_01 = create_dgb("reads/ecoli_sample_reads_01.fasta", k, t)
        print("nb k-mer stored in the graph : ", end='')
        stored1[t] = len(k_mers_ecoli_sample_reads_01)
        print(stored1[t])

        nb_processed = 0
        for key in k_mers_ecoli_sample_reads_01:
            nb_processed += k_mers_ecoli_sample_reads_01[key]
        print("nb k-mer processed is more than : ", end='')
        processed1[t] = nb_processed
        print(nb_processed)


### ecoli_sample_reads ###
print("\n\necoli_sample_reads : ")
stored2 = np.zeros(t_max)
processed2 = np.zeros(t_max)
for t in range(t_max):
    print("\nt = ", end='')
    print(t)

    for k in range(k_min,k_max):
        print("k = ", end='')
        print(k)

        k_mers_ecoli_sample_reads = create_dgb("reads/ecoli_sample_reads.fasta", k, t)
        print("nb k-mer stored in the graph : ", end='')
        stored2[t] = len(k_mers_ecoli_sample_reads)
        print(stored2[t])

        nb_processed = 0
        for key in k_mers_ecoli_sample_reads:
            nb_processed += k_mers_ecoli_sample_reads[key]
        print("nb k-mer processed is more than : ", end='')
        processed2[t] = nb_processed
        print(nb_processed)

### Q5 : plot the frequency ###
t = np.linspace(0, t_max-1, t_max)
fig, ax = plt.subplots(1, 2)
frequency1 = stored1/processed1
frequency2 = stored2/processed2
ax[0].bar(t, frequency1, label="c")
ax[1].bar(t, frequency2, label="d")
ax[0].legend()
ax[1].legend()

plt.show()

