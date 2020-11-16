#!/usr/bin/env python3.6
"""
Find all the apropriate gene families for duplciation dating analysis.

This program will need to be modified for trees not in my chelicerate
project.
"""

import ete3
import glob
import itertools
import re
import sys
import my_module as mod

def get_args():
    """Read user arguments, else prompt the user."""
    if len(sys.argv) == 6:
        return sys.argv[1:]
    elif len(sys.argv) == 5:

        return sys.argv[1:] + ["paralogs_colapsed"]
    else:
        print("\nUSAGE: python3 duplication_species_tree" +
              "gene_tree_dir seq_files_dir outdir *extention\n" +
              "[extention default = paralogs_colapsed]\n")
        exit()

def remove_extensions(tree):
    """Modify the tip labels of a tree so species names remain."""
    for leaf in tree:
        leaf.name = leaf.name[:4]
    return tree

def get_ingroup(tree, outgroups):
    """Get the topology of the ingroup given the list of outgroups."""
    ingroups = []
    for leaf in tree:
        if leaf.name not in outgroups:
            ingroups.append(leaf)
    ingroup_tree = tree.get_common_ancestor(ingroups)
    return ingroup_tree

def add_unique_ids(tree, ingroup_taxa):
    """Add 1 or 2 to ingroup taxa in tree traversal order."""
    counts = {}
    for taxon in ingroup_taxa:
        counts[taxon] = 1
    for leaf in tree:
        if leaf.name in counts:
            name = leaf.name
            leaf.name = name + str(counts[name])
            counts[name] += 1

def correct_n_taxa(tree, ingroup, outgroup):
    """Return true if the number of each taxa fit the specification."""
    if len(tree) < len(ingroup) or len(tree) > len(ingroup) + len(outgroup):
        return 0
    counts = {}
    for taxon in ingroup:
        counts[taxon] = 0
    for taxon in outgroup:
        counts[taxon] = 0
    for leaf in tree:
        counts[leaf.name] += 1
    for taxon in ingroup:
        if counts[taxon] != 2:
            return 0
    outgroups_total = 0
    for taxon in outgroup:
        if counts[taxon] > 1:
            return 0
        outgroups_total += counts[taxon]
    if outgroups_total >= 1:
        return 1
    return 0

def modified_sequence_dict(old_tree, new_tree, seqs, outgroup):
    """Return a seqeunces dict with modified sequence names."""
    #modify seq names to match the treefile
    seq_names = list(seqs.keys())
    for seq in seq_names:
        repl = re.sub("\[|\]", "_", seq)
        seqs[repl] = seqs.pop(seq)

    counts = {}
    
    old_leaves = []
    for leaf in old_tree:
        old_leaves.append(leaf.name)

    new_leaves = []
    for leaf in new_tree:
        new_leaves.append(leaf.name)

    new_seqs = {}
    for i in range(0, len(old_leaves)):
        #get the seqeunce for each leaf and assign it to the
        #right sequence name according to the modified gene tree.
        pattern = re.compile(re.escape(old_leaves[i][5:]))
        key = list(filter(pattern.match, seqs.keys()))[0]
        new_seqs[new_leaves[i]] = seqs[key]

    for taxon in outgroup:
        if taxon not in new_seqs:
            new_seqs[taxon] = "-"

    return new_seqs


def main():
    """
    load a species tree that looks like what I need to emulate
    For each gene family:
        load the tree seqs
        modify the tip labels so that they are just the species names
        check if topologies match.
        if not, check if they match with the deletion of one or two outgroups
        if either of these are true, add them to some sort of list
        write the seqs file but with gene names modified to be just the
        species or species with a number for the duplcited genes.
            This bit might actuall be a bit of pain because i need to ensure
            the 1s are all on the same side, which is the opposite of the 2s
    """
    species_tree, gene_trees, seq_files, outdir, extension = get_args()    
    gene_trees = glob.glob(gene_trees + "/*" + extension)
    
    duplication_tree = ete3.Tree(species_tree)
    pure_sp_tree = remove_extensions(duplication_tree.copy())
    
    outgroup = ["Nstr", "Isca", "Ssca"]
    ingroup = get_ingroup(pure_sp_tree, outgroup)
    ingroup_names = []
    for leaf in ingroup:
        ingroup_names.append(leaf.name)
    ingroup = ingroup_names

    for treefile in gene_trees:
        tree = ete3.Tree(treefile)
        pure_gene_tree = remove_extensions(tree.copy())
        
        #check the number of taxa is to specification
        if not correct_n_taxa(pure_gene_tree, ingroup, outgroup):
            continue

        #add unique names, but only if ingroup is monophyletic
        if pure_gene_tree.check_monophyly(values=ingroup_names, target_attr="name"):
            add_unique_ids(pure_gene_tree, ingroup_names)
        else: #here the ingroup is not monophyletic so skip
            continue
        
        #Continue of rf distance is 0. Given our previous taxa occupency 
        #checks this well confirm our trees are correct
        if pure_gene_tree.robinson_foulds(duplication_tree)[0] == 0:
            #find sequence file
            base_name = treefile.split("/")[1].split("_")[0]
            seqs = mod.read_fasta(seq_files + "/" + base_name + ".fa")
            
            new_seqs = modified_sequence_dict(tree, pure_gene_tree, seqs, outgroup)
            
            out = open(outdir + "/" + base_name + ".fa", "w")
            for key, value in new_seqs.items():
                out.write(">" + key + "\n" + value + "\n")
            out.close()
        

if __name__ == "__main__":
    main()
