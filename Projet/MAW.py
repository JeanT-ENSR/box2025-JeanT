#! /usr/bin/env python3

from xopen import xopen
from readfa import readfq
import numpy as np
import sys
from collections import defaultdict

def open(fastaFilePath):
    sequences = []
    with xopen(fastaFilePath) as fa:
        for _,seq,_ in readfq(fa):
            sequences.append(seqfa)
    return sequences

