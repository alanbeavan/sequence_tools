#!usr/bin/env python3.6
"""Take a random selection of sites and make a new superalignment file."""

import sys
import random
import my_module as mod

def get_args():
    """Get user arguments."""
    if len(sys.argv) == 4:
        return sys.argv[1:]
    print("USAGE: python3 jackknife.py infile, outfile, n_sites\n")
    exit()

def main():
    """Do the things."""
    infile, outfile, n_sites = get_args()
    seqs = mod.read_fasta(infile)
    flag = 1
    out = open(outfile, "w")
    for key, value in seqs.items():
        if flag:
            indicies = sorted(random.sample(range(len(value)), int(n_sites)))
            flag = 0
        out.write(">" + key + "\n" + "".join([value[i] for i in indicies]) + "\n")
    out.close()

if __name__ == "__main__":
    main()
