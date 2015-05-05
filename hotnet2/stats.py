# -*- coding: iso-8859-1 -*-
from collections import defaultdict
import multiprocessing as mp
import networkx as nx
import hotnet2 as hn
import hnio

strong_ccs = nx.strongly_connected_components

def num_components_min_size(G, sizes):
    """Return a list of the number of connected components of size at least s for each s in sizes. 
    
    Arguments:
    G -- a networkx Graph or DiGraph
    sizes -- an iterable of minimum connected component sizes
    
    """
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    cc_sizes = [len(cc) for cc in ccs]
    return [sum(1 for cc_size in cc_sizes if cc_size >= s) for s in sizes]

def significance_wrapper((infmat, index2gene, heat_permutation, delta, sizes, directed)):
    sim, index2gene = hn.similarity_matrix(infmat, index2gene, heat_permutation, directed)
    G = hn.weighted_graph(sim, index2gene, delta, directed)
    return num_components_min_size(G, sizes)

def network_significance_wrapper((network_path, infmat_name, index2gene, heat, delta, sizes, directed)):
    permuted_mat = hnio.load_hdf5(network_path)[infmat_name]
    sim, index2gene = hn.similarity_matrix(permuted_mat, index2gene, heat, directed)
    G = hn.weighted_graph(sim, index2gene, delta, directed)
    return num_components_min_size(G, sizes)

def calculate_permuted_cc_counts_network(network_paths, infmat_name, index2gene, heat, delta,
                                         sizes=range(2,11), directed=True, num_cores=1):
    """Return a dict mapping a CC size to a list of the number of CCs of that size or greater in
    each permutation.
    """
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map
    
    args = [(network_path, infmat_name, index2gene, heat, delta, sizes, directed)
            for network_path in network_paths] 
    all_counts = map_fn(network_significance_wrapper, args)
    
    if num_cores != 1:
        pool.close()
        pool.join()

    # Parse the results into a map of k -> counts
    size2counts = defaultdict(list)
    for counts in all_counts:
        for size, count in zip(sizes, counts): size2counts[size].append(count)            

    return size2counts

def calculate_permuted_cc_counts(infmat, index2gene, heat_permutations, delta,
                                 sizes=range(2,11), directed=True, num_cores=1):
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
    directed -- whether the graph constructed from each permuted similarity matrix should be directed
    num_cores -- number of cores to use for running in parallel
    
    """
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map
    
    args = [(infmat, index2gene, heat_permutation, delta, sizes, directed)
            for heat_permutation in heat_permutations] 
    all_counts = map_fn(significance_wrapper, args)
    
    if num_cores != 1:
        pool.close()
        pool.join()

    # Parse the results into a map of k -> counts
    size2counts = defaultdict(list)
    for counts in all_counts:
        for size, count in zip(sizes, counts): size2counts[size].append(count)            

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
    size2stats = dict()
    for size, counts in size2counts_permuted.items():
        observed = size2counts_real[size]
        expected = sum(counts) / num_permutations
        pval = len([c for c in counts if c >= observed]) / num_permutations
        size2stats[size] = dict(observed=observed, expected=expected, pval=pval)
    
    return size2stats
