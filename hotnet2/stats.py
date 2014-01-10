# -*- coding: iso-8859-1 -*-
import hotnet2 as hn
import networkx as nx
import multiprocessing as mp
import scipy.io
strong_ccs = nx.strongly_connected_components

def num_components_min_size(G, sizes):
    """Return a list of the number of connected components of size at least s for each s in sizes. 
    
    Arguments:
    G -- a networkx Graph or DiGraph
    sizes -- an iterable of minimum connected component sizes
    
    """
    return [len([ cc for cc in strong_ccs(G) if len(cc) >= s]) for s in sizes]

def significance_wrapper((infmat, index2gene, heat_permutation, delta, sizes)):
    M, index2gene = hn.induce_infmat(infmat, index2gene, sorted(heat_permutation.keys()))
    h = hn.heat_vec(heat_permutation, index2gene)
    sim_mat = hn.similarity_matrix(M, h)
    G = hn.weighted_graph(sim_mat, index2gene, delta)
    return num_components_min_size(G, sizes)

def calculate_permuted_cc_counts(infmat, index2gene, heat_permutations, delta,
                                 sizes=range(2,11), parallel=True):
    """Return a dict mapping a CC size to a list of the number of CCs of that size or greater in
    each permutation.
    
    Arguments:
    infmat -- 2D ndarray representing an influence matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    heat_permutations -- iterable of dicts mapping a gene name to the permuted heat score for that gene
    delta -- threshold for edge weight removal
    sizes -- list of sizes for which the number of connected components of that sizes should be
             calculated in each permutation
    parallel -- whether the similarity matrix and connected components for each permutation should
                be calculated in parallel
    
    """
    if parallel:
        pool = mp.Pool(25)
        map_fn = pool.map
    else:
        map_fn = map
    
    args = [(infmat, index2gene, heat_permutation, delta, sizes)
            for heat_permutation in heat_permutations] 
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
    """Return a dict mapping a CC size to a tuple with the expected number of CCs of at least that
    size based on permuted data, the observed number of CCs of at least that size in the real data,
    and the p-value for the observed number.
    
    sizes2counts_real -- dict mapping a CC size to the number of CCs of that size or greater
                         observed in the real data
    sizes2counts_permuted -- dict mapping a CC size to a list of the number of CCs of that size or
                             greater in each permutation
    num_permutations -- the number of permuted data sets represented in sizes2counts_permuted
     
    """
    num_permutations = float(num_permutations)
    size2stats = dict([(s, []) for s in size2counts_permuted.keys()])
    for size, counts in size2counts_permuted.items():
        observed = size2counts_real[size]
        expected = sum(counts) / num_permutations
        pval = len([ c for c in counts if c >= observed ]) / num_permutations
        size2stats[size] = dict(observed=observed, expected=expected, pval=pval)
    
    return size2stats
