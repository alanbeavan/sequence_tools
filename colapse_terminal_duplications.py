#!/usr/bin/env python3.6
"""Remove one of a pair of inparalogs from a gene tree.


This is written for orhtofinder output. Modify it for your trees."""

import ete3
import re
import sys
import my_module as mod

def get_args():
    """Get user arguments or prompt the user."""
    if len(sys.argv) == 3:
        return sys.argv[1:]
    else:
        print("USAGE python3 colapse_terminal_duplications.py treefile seqfile")
        exit()

def main():
    """
    Remove the shortest paralog from a terminal duplicaiton.
    
    Take a tree and seqeunce file from user arguments.
    While there are terminal duplicaitons (2 seqeunces from the same species
    who are eachothers nearest relative):
        Remove the sequence of the tuple with the shortest sequence.
    """
    treefile, seqfile = get_args()
    seqs = mod.read_fasta(seqfile)

    #modify seq names to match the treefile
    seq_names = list(seqs.keys())
    for seq in seq_names:
        repl = re.sub("\[|\]", "_", seq)
        seqs[repl] = seqs.pop(seq)


    tree = ete3.Tree(treefile, format = 1)
    #for node in nodes, if both descendants are a tip from the same species, delete the right one and start again
    done = 0
    while done == 0 and len(tree) > 1:
        print(tree)
        old_tree = tree.copy()
        for node in tree.traverse():
            descendants = list(node.iter_descendants())
            if len(node) == 1 and not node.is_leaf():
                node.delete()
                break
            if len(node) == 2:
                if descendants[0].name[:4] == descendants[1].name[:4]:

                    names = []
                    lengths = {}
                    for descendant in descendants:
                        names.append(descendant.name)
                        lengths[descendant.name] = len(seqs[descendant.name[5:]])
                    
                    to_delete = min(names, key=lengths.__getitem__)
                    to_delete = tree.search_nodes(name = to_delete)[0]
                    to_delete.detach()
                    break
                    done = 1
        if len(list(old_tree.traverse())) == len(list(tree.traverse())):
            done = 1

    #Write tree
    if len(tree) > 1:
        tree.write(outfile = treefile + ".paralogs_colapsed")
        out = open(seqfile + ".paralogs_colapsed", "w")
        for leaf in tree:
            pattern = re.compile(leaf.name[5:])
            key = list(filter(pattern.match, seqs.keys()))[0]
            out.write(">" + leaf.name[5:] + "\n" + seqs[key])
        out.close()


if __name__ == "__main__":
    main()
