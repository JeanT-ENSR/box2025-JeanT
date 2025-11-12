import numpy as np
# import mmh3 # librairie de fonctions de hash

nb_hashes = 7
size_max = 100000000

item = "ACGGACGACGACT"
#for seed in range(nb_hashes):
#    key = mmh3.hash(item, seed, signed=false) % size_max
#    print(f"Seed {seed} is {key}")

def hash(item, seed):
    sum = 0
    for i in range(len(item)):
        v = ord(item[i])
        sum = (sum + v)*seed
    return sum

bloom_array = np.zeros(size_max, dtype=np.int8)

def insert_counting(item):
    for seed in range(nb_hashes):
        key = hash(item,seed) % size_max
        bloom_array[key] += 1
    return

def query_counting(item):
    key = hash(item,0) % size_max
    min_val = bloom_array[key]
    for seed in range(1,nb_hashes):
        key = hash(item,seed) % size_max
        min_val = min(min_val,bloom_array[key])
    return min_val

insert_counting(item)
print(query_counting(item))
print(query_counting("ACAC"))