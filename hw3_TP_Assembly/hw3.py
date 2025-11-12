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
k_min = 31
k_max = 32
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
frequency1 = processed1/stored1
frequency2 = processed2/stored2
ax[0].bar(t, frequency1, label="c")
ax[1].bar(t, frequency2, label="d")
ax[0].legend()
ax[1].legend()

plt.show()


### Unitigs ###
def unitig_from(dgb, k, s):
    unitig = s
    if not (s in dgb):
        return unitig
    else:
        # On avance dans le graph à partir de s
        continu = True
        nb = 1
        while continu :
            possibility = 0
            a = (unitig[nb:(30+nb)] + "A") in dgb
            c = (unitig[nb:(30+nb)] + "C") in dgb
            g = (unitig[nb:(30+nb)] + "G") in dgb
            t = (unitig[nb:(30+nb)] + "T") in dgb
            if a :
                possibility = possibility + 1
            if c :
                possibility = possibility + 1
            if g :
                possibility = possibility + 1
            if t :
                possibility = possibility + 1

            if possibility == 1:
                nb = nb + 1
                if a :
                    unitig = unitig + "A"
                elif c :
                    unitig = unitig + "C"
                elif g :
                    unitig = unitig + "G"
                elif t :
                    unitig = unitig + "T"
            else:
                continu = False

        # On recule dans le graph à partir de s
        continu = True
        while continu :
            possibility = 0
            a = ("A" + unitig[0:30]) in dgb
            c = ("C" + unitig[0:30]) in dgb
            g = ("G" + unitig[0:30]) in dgb
            t = ("T" + unitig[0:30]) in dgb
            if a :
                possibility = possibility + 1
            if c :
                possibility = possibility + 1
            if g :
                possibility = possibility + 1
            if t :
                possibility = possibility + 1

            if possibility == 1:
                if a :
                    unitig = "A" + unitig
                elif c :
                    unitig = "C" + unitig
                elif g :
                    unitig = "G" + unitig
                elif t :
                    unitig = "T" + unitig
            else:
                continu = False

        return unitig


k = 31
size_unitig = np.zeros((4,5))
for t in range(1,5):
    print("\nt = ", end='')
    print(t)

    print("ecoli_genome_150k.fa :")
    unitig = unitig_from( (create_dgb("ecoli_genome_150k.fa",k,t)), k, 'CGCTCTGTGTGACAAGCCGGAAACCGCCCAG')
    print("len( ", end='')
    print(unitig, end='')
    print(") = ", end='')
    size_unitig[t-1][0] = len(unitig)
    print(size_unitig[t-1][0])

    print("\nreads/ecoli_sample_perfect_reads_forward.fasta :")
    unitig = unitig_from((create_dgb("reads/ecoli_sample_perfect_reads_forward.fasta",k,t)), k, 'CGCTCTGTGTGACAAGCCGGAAACCGCCCAG')
    print("len( ", end='')
    print(unitig, end='')
    print(") = ", end='')
    size_unitig[t-1][1] = len(unitig)
    print(size_unitig[t-1][1])

    print("\nreads/ecoli_sample_perfect_reads.fasta :")
    unitig = unitig_from((create_dgb("reads/ecoli_sample_perfect_reads.fasta",k,t)), k, 'CGCTCTGTGTGACAAGCCGGAAACCGCCCAG')
    print("len( ", end='')
    print(unitig, end='')
    print(") = ", end='')
    size_unitig[t-1][2] = len(unitig)
    print(size_unitig[t-1][2])

    print("\nreads/ecoli_sample_reads_01.fasta :")
    unitig = unitig_from((create_dgb("reads/ecoli_sample_reads_01.fasta",k,t)), k, 'CGCTCTGTGTGACAAGCCGGAAACCGCCCAG')
    print("len( ", end='')
    print(unitig, end='')
    print(") = ", end='')
    size_unitig[t-1][3] = len(unitig)
    print(size_unitig[t-1][3])

    print("\nreads/ecoli_sample_reads.fasta :")
    unitig = unitig_from((create_dgb("reads/ecoli_sample_reads.fasta",k,t)), k, 'CGCTCTGTGTGACAAGCCGGAAACCGCCCAG')
    print("len( ", end='')
    print(unitig, end='')
    print(") = ", end='')
    size_unitig[t-1][4] = len(unitig)
    print(size_unitig[t-1][4])

print("\nt = 1 :")
print(size_unitig[0][0], end=' ')
print(size_unitig[0][1], end=' ')
print(size_unitig[0][2], end=' ')
print(size_unitig[0][3], end=' ')
print(size_unitig[0][4])

print("\nt = 2 :")
print(size_unitig[1][0], end=' ')
print(size_unitig[1][1], end=' ')
print(size_unitig[1][2], end=' ')
print(size_unitig[1][3], end=' ')
print(size_unitig[1][4])

print("\nt = 3 :")
print(size_unitig[2][0], end=' ')
print(size_unitig[2][1], end=' ')
print(size_unitig[2][2], end=' ')
print(size_unitig[2][3], end=' ')
print(size_unitig[2][4])

print("\nt = 4 :")
print(size_unitig[3][0], end=' ')
print(size_unitig[3][1], end=' ')
print(size_unitig[3][2], end=' ')
print(size_unitig[3][3], end=' ')
print(size_unitig[3][4])

