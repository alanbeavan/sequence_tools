#!/usr/bin/env python3.6

import sys, os, subprocess, re
sys.path.append("/Users/ab17362/OneDrive - University of Bristol/Python_modules")
from my_module import *

#takes a fasta as input and loops through the blast databases in a directory, checking if the best hit is for each gene in the gene family
#corresponds to the representative of that species that is present in the fasta file

if len(sys.argv) != 3:
    print("Usage:  python3 fasta_file db_directory\n\ndb_directory can be full or relative path. We're easy")
    exit()
else:
    seq_file = sys.argv[1]
    db_directory = sys.argv[2]

seqs = get_file_data(seq_file)
files = os.listdir(db_directory)
dbs = sorted(set([i.split(".", 1)[0] for i in files]))

for db in dbs:
    subprocess.run("blastp -query " + seq_file + " -db " + db_directory + "/"  + db + " -out out_temp -outfmt 6 -evalue 1e-09 -num_alignments 1", shell = True)
    lines = get_file_data("out_temp")
    for line in lines:
        hit = line.split("\t")[1]
        r = re.compile(hit)
        if len(list(s for s in seqs if hit in s)) >= 1:
#            print(line)
            next
        else:
#            print(line)
            print(seq_file + " reciprocicity not fullfilled")
            os.remove("out_temp")
            exit()
    os.remove("out_temp")
print(seq_file + " reciprocicity fullfilled")
