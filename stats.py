#!/usr/bin/python

# Import required modules
import numpy as np, sys
from hotnet2 import *
from networkx import strongly_connected_components as strong_ccs

def num_components_min_size(G, sizes):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    return [len([ cc for cc in ccs if len(cc) > s]) for s in sizes]


def significance_wrapper( (infmat, infmat_index, genes, h, delta, filtered_indices, sizes, directed, i) ):
    # print "\t\t- Permutation", i
    permuted_genes = permute_and_filter_genes(genes, filtered_indices)
    M, gene_index, _ = induce_infmat( infmat, infmat_index, permuted_genes)
    sim_mat, _ = similarity_matrix(M, h, gene_index, directed)
    G = weighted_graph(sim_mat, gene_index, delta)
    return num_components_min_size( G, sizes )


def permute_and_filter_genes( genes, filtered_indices ):
    permuted_genes = list(genes)
    shuffle( permuted_genes )
    return [ permuted_genes[i] for i in filtered_indices ]


from random import shuffle
import multiprocessing as mp
def calculate_permuted_cc_counts(infmat, infmat_index, genes, h, delta, filtered_genes, n,
                                 sizes=range(2, 11), directed=True, parallel=True):
    # Find indices of filtered genes in the set of all genes
    filtered_indices = [genes.index(g) for g in filtered_genes]
    
    # Report parameters of run
    print "* Performing permuted heat statistical signifcance..."
    print "\t- Using no. of components >= k (k \\in",
    print "[%s, %s]) as statistic" % (min(sizes), max(sizes))
    print "\t- Running permutations:"
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map
        
    args = [(infmat, infmat_index, genes, h, delta, filtered_indices, sizes, directed, i) for i in range(n)]
    all_counts = map_fn(significance_wrapper, args)
    
    if parallel:
        pool.close()
        pool.join()

    # Parse the results into a map of k -> counts
    size2counts = dict([(s, []) for s in sizes])
    for counts in all_counts:
        for size, count in zip(sizes, counts): size2counts[size].append( count )            

    return size2counts


def compute_statistics(size2counts_real, size2counts_permuted, num_permutations):
    size2stats = dict([(s, []) for s in size2counts_permuted.keys()])
    for size, counts in size2counts_permuted.items():
        observed = size2counts_real[size]
        expected = sum(counts) / num_permutations
        pval = len([ c for c in counts if c >= observed ]) / float( num_permutations )
        size2stats[size] = dict(observed=observed, expected=expected, pval=pval)
    
    return size2stats
