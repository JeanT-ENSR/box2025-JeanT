import numpy as np
import mmh3 # librairie de fonctions de hash

nb_hashes = 7
size_max = 100000000
item = "ACGGACGACGACT"
for seed in range(nb_hashes):
    key = mmh3.hash(item, seed, signed=false) % size_max
    print(f"Seed {seed} is {key}")