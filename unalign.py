#!/usr/bin/env python3.6

import sys, re
sys.path.append("/Users/ab17362/OneDrive - University of Bristol/Python_modules")
from my_module import *

#just removes "-" characters from sequences (fasta)
if len(sys.argv) != 3:
    print("Usage\n\npython3 unalign.py infile outfile")
    exit()


lines = get_file_data(sys.argv[1])
newlines = []
seq = ""
for line in lines:
    if line.startswith(">"):
        if seq:
            newlines.append(seq)
        newlines.append(line)
        seq = ""
    else:
        to_add = re.sub("-", "", line)
        seq = seq + to_add

f = open(sys.argv[2], "w")
f.write("\n".join(newlines))
f.close()
