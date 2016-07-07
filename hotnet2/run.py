#!/usr/bin/env python

# Load required modules
import json
import sys, os, shutil
import numpy as np

# Load local modules
import heat as hnheat, hotnet2 as hn, hnio, stats, permutations as p
from delta import get_deltas_for_network, get_deltas_for_heat
from constants import *

def run_helper(args, infmat, full_index2gene, G, nname, pnp, heat, hname, addtl_genes, get_deltas_fn, infmat_name="PPR", max_cc_sizes=[5, 10, 15, 20], verbose=0):
    """Helper shared by runHotNet2 and runClassicHotNet.
    """
    # Perform delta selection (if necessary)
    if args.deltas:
        deltas = args.deltas
    else:
        deltas = get_deltas_fn(full_index2gene, heat, args.network_permutations, args.num_cores, infmat, addtl_genes, pnp, infmat_name, max_cc_sizes, verbose)

    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, True, verbose=verbose)

    results = []
    for delta in deltas:

        # find connected components
        G = hn.weighted_graph(sim, index2gene, delta, directed=True)
        ccs = hn.connected_components(G, args.min_cc_size)

        # calculate significance (using all genes with heat scores)
        if verbose > 4:
            print "* Performing permuted heat statistical significance..."
            print "\t- Using no. of components >= k (k \\in",
            print "[%s, %s]) as statistic" % (min(HN2_STATS_SIZES), max(HN2_STATS_SIZES))

        heat_permutations = p.permute_heat(heat, full_index2gene.values(),
                                           args.heat_permutations, addtl_genes,
                                           args.num_cores)
        sizes2counts = stats.calculate_permuted_cc_counts(infmat, full_index2gene,
                                                          heat_permutations, delta, HN2_STATS_SIZES, True,
                                                          args.num_cores)
        real_counts = stats.num_components_min_size(G, HN2_STATS_SIZES)
        size2real_counts = dict(zip(HN2_STATS_SIZES, real_counts))
        sizes2stats = stats.compute_statistics(size2real_counts, sizes2counts,
                                               args.heat_permutations)

        # sort ccs list such that genes within components are sorted alphanumerically, and components
        # are sorted first by length, then alphanumerically by name of the first gene in the component
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)

        # Record the results for this delta
        results.append( (ccs, sizes2stats, delta) )

    return results

def get_deltas_hotnet2(full_index2gene, heat, num_perms, num_cores, _infmat, _addtl_genes,
                       permuted_networks_path, infmat_name, max_cc_sizes, verbose):
    # find smallest delta
    deltas = get_deltas_for_network(permuted_networks_path, heat, infmat_name, full_index2gene,
                                       MAX_CC_SIZE, max_cc_sizes, False, num_perms, num_cores, verbose)

    # and run HotNet with the median delta for each size
    return [np.median(deltas[size]) for size in deltas]

def get_deltas_classic(full_index2gene, heat, num_perms, num_cores, infmat, addtl_genes, min_cc_size, max_cc_size, verbose):
    # find delta that maximizes # CCs of size >= min_cc_size for each permuted data set
    deltas = get_deltas_for_heat(infmat, full_index2gene, heat, addtl_genes, num_perms, NUM_CCS,
                                    [min_cc_size], True, num_cores, verbose)

    # find the multiple of the median delta s.t. the size of the largest CC in the real data
    # is <= MAX_CC_SIZE
    medianDelta = np.median(deltas[min_cc_size])

    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, False)

    for i in range(1, 11):
        G = hn.weighted_graph(sim, index2gene, i*medianDelta)
        largest_cc_size = max([len(cc) for cc in hn.connected_components(G)])
        if largest_cc_size <= max_cc_size:
            break

    # and run HotNet with that multiple and the next 4 multiples
    return [i*medianDelta for i in range(i, i+5)]
