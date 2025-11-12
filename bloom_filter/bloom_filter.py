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

bloom_array = np.zeros(size_max, dtype=bool)

def insert(item):
    for seed in range(nb_hashes):
        key = hash(item,seed) % size_max
        bloom_array[key] = True
    return

def query(item):
    is_in = True
    for seed in range(nb_hashes):
        key = hash(item,seed) % size_max
        is_in = is_in and bloom_array[key]
    return is_in

insert(item)
print(query(item))
print(query("ACAC"))