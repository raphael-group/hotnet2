# -*- coding: iso-8859-1 -*-
import hotnet2 as hn
import networkx as nx
import multiprocessing as mp
strong_ccs = nx.strongly_connected_components

def num_components_min_size(G, sizes):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    return [len([ cc for cc in ccs if len(cc) > s]) for s in sizes]

def significance_wrapper((infmat, index2gene, heat_permutation, delta, sizes, directed)):
    M, index2gene = hn.induce_infmat(infmat, index2gene, sorted(heat_permutation.keys()))
    h = hn.heat_vec(heat_permutation, index2gene)
    sim_mat = hn.similarity_matrix(M, h, directed)
    G = hn.weighted_graph(sim_mat, index2gene, delta, directed)
    return num_components_min_size(G, sizes)

def calculate_permuted_cc_counts(infmat, index2gene, heat_permutations, delta,
                                 sizes=range(2,11), directed=True, parallel=True):
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
    
    args = [(infmat, index2gene, heat_permutation, delta, sizes, directed)
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
    size2stats = dict([(s, []) for s in size2counts_permuted.keys()])
    for size, counts in size2counts_permuted.items():
        observed = size2counts_real[size]
        expected = sum(counts) / num_permutations
        pval = len([ c for c in counts if c >= observed ]) / float( num_permutations )
        size2stats[size] = dict(observed=observed, expected=expected, pval=pval)
    
    return size2stats
