#!/usr/bin/env python3.6

import sys, glob
sys.path.append("/Users/ab17362/OneDrive - University of Bristol/Python_modules")
from my_module import *

species = []
species_lines = get_file_data("species_file")
for line in species_lines:
    species.append(line)


seq_files = glob.glob("SEQ_FILES")
for file in seq_files:
    seqs = read_fasta(file)
    f = open("_".join(file.split("_")[0:2]) + "_complete.fa", "w")
    for s in species:
        if s in seqs.keys():
            f.write(">" + s + "\n" + seqs[s] + "\n")
        else:
            f.write(">" + s + "\n-\n")
