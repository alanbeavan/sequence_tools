#!/usr/bin/env python3.6
"""Extract the coding sequences from a transcriptome."""

import sys
import my_module as mod

def get_args():
    """Get user arguments for infile and outfile."""
    if len(sys.argv) == 3:
        return sys.argv[1:]
    print("USAGE:  python3 extract_cds.py infile outfile\n\n")
    exit()

def remove_guff(seqs):
    """Remove anything before the first start codon and after stop codon."""
    new_seqs = {}
    stop_codons = ["TGA", "TAA", "TAG"]
    for key, value in seqs.items():
        new_seq = ""
        for i in range(len(value)-2):
            if value[i:i+3] == "ATG":
                break

        for j in range(i, len(value)-2, 3):
            if value[j:j+3] in stop_codons:
                new_seqs[key] = value[i:j+3]
                break

    return new_seqs



def main():
    """
    Perform the following.
    
    1. Get user input of infile and outfile
    2. Read in the sequences from the transcriptome
    3. For each sequence, take the first start codon in any frame and read
    through until a stop codon is reached. Write this to the outfile. If a
    stop codon is not reached, just remove the sequence from the outfile.
    """
    infile, outfile = get_args()
    seqs = mod.read_fasta(infile)

    coding_seqs = remove_guff(seqs)
    
    out = open(outfile, "w")
    for key, value in coding_seqs.items():
        out.write(">" + key + "\n" + value + "\n")
    out.close()

if __name__ == "__main__":
    main()
