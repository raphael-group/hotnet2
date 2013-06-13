#!/usr/bin/python

# Import required modules
import numpy as np, sys
from data import load_index
from hotnet2 import *
from networkx import strongly_connected_components as strong_ccs

def num_components_min_size( DiG, sizes ):
    ccs = strong_ccs( DiG )
    return [len([ cc for cc in ccs if len(cc) > s]) for s in sizes]

def significance_wrapper( (infmat, infmat_index, genes, h, delta, filtered_indices, sizes, i) ):
    # print "\t\t- Permutation", i
    permuted_genes = permute_and_filter_genes( genes, filtered_indices )
    M, gene_index, _ = induce_infmat( infmat, infmat_index, permuted_genes )
    sim_mat, _ = similarity_matrix( M, h, gene_index )
    DiG = weighted_graph( sim_mat, gene_index, delta )
    return num_components_min_size( DiG, sizes )

def permute_and_filter_genes( genes, filtered_indices ):
    permuted_genes = list(genes)
    shuffle( permuted_genes )
    return [ permuted_genes[i] for i in filtered_indices ]

from random import shuffle
import multiprocessing as mp
def compute_significance( DiG, infmat, infmat_index, genes, h, delta,
                          filtered_genes, n, sizes=range(2, 11), parallel=True):
    # Find indices of filtered genes in the set of all genes
    filtered_indices = [genes.index(g) for g in filtered_genes]
    
    # Report parameters of run
    print "* Performing permuted heat statistical signifcance..."
    print "\t- Using no. of components >= k (k \\in",
    print "[%s, %s]) as statistic" % (min(sizes), max(sizes))
    print "\t- Running permutations:"
    if parallel:
        # Submit jobs to all available processors
        pool = mp.Pool()
        args = [(infmat, infmat_index, genes, h, delta, filtered_indices, sizes, i )
                for i in range(n)]
        all_counts = pool.map( significance_wrapper, args )
        pool.close()
        pool.join()

    else:
        all_counts = []
        for i in range(1, n+1):
            print "\t\t- Permutation", i
            permuted_genes = permute_and_filter_genes( genes, filtered_indices )
            M, gene_index, _ = induce_infmat( infmat, infmat_index, permuted_genes )
            sim_mat, _ = similarity_matrix( M, h, gene_index )
            permuted_DiG = weighted_graph( sim_mat, gene_index, delta )
            counts = num_components_min_size( permuted_DiG, sizes )
            all_counts.append( counts )

    # Parse the results into a map of k -> counts
    size2counts = dict([(s, []) for s in sizes])
    for counts in all_counts:
        for size, count in zip(sizes, counts): size2counts[size].append( count )            

    # Compute the scores on the given network
    print "\t- Counts for given delta on real data:"
    real_counts = num_components_min_size( DiG, sizes )
    size2real_counts = dict(zip(sizes, real_counts))
    for size, count in size2real_counts.items():
        print "\t\tk=%s -> %s" % (size, count)

    # Compute statistics from size2stats
    size2stats = dict([(s, []) for s in sizes])
    for size, counts in size2counts.items():
        observed = size2real_counts[size]
        expected = sum(counts) / n
        pval = len([ c for c in counts if c >= observed ]) / float( n )
        size2stats[size] = dict(observed=observed, expected=expected, pval=pval)
    
    return size2stats
