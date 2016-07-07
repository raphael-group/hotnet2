#!/usr/bin/env python

import sys, os, networkx as nx, numpy
from itertools import product, combinations
from collections import defaultdict
from run import run_helper, get_deltas_hotnet2
from constants import *
from heat import filter_heat, filter_heat_to_network_genes
from permutations import permute_heat

def count_consensus(consensus, sizes=HN2_STATS_SIZES):
    cc_sizes = [ len(d['core'] + d['expansion']) for d in consensus ]
    return dict( (s, sum(1 for cc_size in cc_sizes if cc_size >= s)) for s in sizes )

def consensus_with_stats(args, networks, heats, verbose=0):
    # Run with the input heat
    single_runs, consensus, linkers, auto_deltas = consensus_run( args, networks, heats, verbose )

    # Generate permuted heats
    np = args.consensus_permutations
    permuted_single_runs = defaultdict(list)
    for (infmat, indexToGene, G, nname, pnp), (heat, hname) in product(networks, heats):
        # 1) Filter the heat scores
        # 1a) Remove enes not in the network
        heat = filter_heat_to_network_genes(heat, set(indexToGene.values()), verbose)

        # 1b) Genes with score 0 cannot be in output components, but are eligible for heat in permutations
        heat, addtl_genes = filter_heat(heat, None, False, 'There are ## genes with heat score 0')

        for permutation in permute_heat(heat, indexToGene.values(), np, addtl_genes, args.num_cores):
            result = run_helper(args, infmat, indexToGene, G, nname, pnp, heat, hname, addtl_genes, get_deltas_hotnet2, HN2_INFMAT_NAME, HN2_MAX_CC_SIZES, verbose=verbose)
            permuted_single_runs[(hname, nname)].append(result)

    # Run consensus to compute observed statistics
    network_heat_pairs = permuted_single_runs.keys()
    permuted_counts = []
    for i in range(args.heat_permutations):
        runs = [ (n, h, permuted_single_runs[(n, h)][i]) for n, h in network_heat_pairs ]
        permuted_consensus, _, _ = identify_consensus( runs, verbose=verbose )
        permuted_counts.append(count_consensus(permuted_consensus))

    # Summarize stats
    consensus_stats = dict()
    for k, count in count_consensus(consensus).iteritems():
        empirical = [ permuted_count[k] for permuted_count in permuted_counts ]
        if np == 0:
            pval     = 1.
            expected = 0.
        else:
            expected = numpy.mean(empirical)
            pval     = sum(1. for p in empirical if p >= count )/np
        consensus_stats[k] = dict(observed=count, expected=expected, pval=pval)

    return single_runs, consensus, linkers, auto_deltas, consensus_stats

def consensus_run(args, networks, heats, verbose):
    # Perform the single runs
    single_runs = []
    for (infmat, indexToGene, G, nname, pnp), (heat, hname) in product(networks, heats):
        # Simple progress bar
        if args.verbose > 0: print '\t-', nname, hname

        # 1) Filter the heat scores
        # 1a) Remove enes not in the network
        heat = filter_heat_to_network_genes(heat, set(indexToGene.values()), verbose)

        # 1b) Genes with score 0 cannot be in output components, but are eligible for heat in permutations
        heat, addtl_genes = filter_heat(heat, None, False, 'There are ## genes with heat score 0')

        if args.verbose > 1:
            print "\t\t- Loaded '%s' heat scores for %s genes" % (hname, len(heat))

        result = run_helper(args, infmat, indexToGene, G, nname, pnp, heat, hname, addtl_genes, get_deltas_hotnet2, HN2_INFMAT_NAME, HN2_MAX_CC_SIZES, args.verbose)
        single_runs.append( (nname, hname, result) )

    # Perform the consensus
    consensus, linkers, auto_deltas = identify_consensus( single_runs, verbose=verbose )

    return single_runs, consensus, linkers, auto_deltas

def identify_consensus(single_runs, pval_threshold=0.01, min_cc_size=2, verbose=0):
    if verbose > 0: print '* Constructing HotNet(2) consensus network...'

    # Choose single runs and count the number of networks
    components, networks, auto_deltas = choose_deltas(single_runs, pval_threshold, min_cc_size, verbose)
    num_networks = len(set(networks))

    # Create the consensus graph
    edges = consensus_edges(components, networks)
    G = nx.Graph()
    G.add_weighted_edges_from( (u, v, w) for (u, v), w in edges.iteritems() )

    # Extract the connected components when restricted to edges in all networks.
    H = nx.Graph()
    H.add_edges_from( (u, v) for (u, v), w in edges.iteritems() if w >= num_networks )
    consensus = [ set(cc) for cc in nx.connected_components( H ) ]
    consensus_genes = set( g for cc in consensus for g in cc )

    # Expand each consensus by adding back any edges not in all networks.
    expanded_consensus = []
    linkers = set()
    for cc in consensus:
        other_consensus_genes = consensus_genes - cc
        neighbors = set( v for u in cc for v in G.neighbors(u) if v not in consensus_genes )
        expansion = set()
        for u in neighbors:
            cc_networks = max( G[u][v]['weight'] for v in set(G.neighbors(u)) & cc )
            consensus_neighbors = set( v for v in G.neighbors(u) if v in other_consensus_genes and G[u][v]['weight'] >= cc_networks )
            if len(consensus_neighbors) > 0:
                if any([ G[u][v]['weight'] > 1 for v in consensus_neighbors ]):
                    linkers.add( u )
            else:
                expansion.add( u )
        expanded_consensus.append( dict(core=list(cc), expansion=list(expansion)) )

    consensus_genes = set( g for cc in expanded_consensus for g in cc['core'] + cc['expansion'] )
    linkers -= consensus_genes

    return expanded_consensus, linkers, auto_deltas

# Choose the results files when given directories of results files.
def choose_deltas(single_runs, pval_threshold, min_cc_size, verbose=0):
    components, networks, auto_deltas = [], [], []
    # Consider each run.
    for network, heat, run in single_runs:
        result_statistics = []
        for ccs, stats, delta in run:
            # Consider each \delta value for each run.
            # For each \delt value, load the results for each \delta value and record the
            # number of component sizes k with p-values below a given threshold when k \geq
            # min_cc_size.
            count = sum( 1 for k in stats if int(k) >= min_cc_size and stats[k]['pval'] > pval_threshold )
            result_statistics.append((count, delta, ccs))

        # Find the smallest \delta value with the largest number of component sizes k with p-values
        # below a threshold. For convenience, we find the smallest number of \delta values greater
        # greater than or equal to the threshold, and we sort from smallest to largest.
        selected_result_statistics = sorted(result_statistics)[0]
        if verbose > 3:
            print '\t- Selected delta = {} for {}...'.format(selected_result_statistics[1], directory)

        components.append(selected_result_statistics[2])
        auto_deltas.append(selected_result_statistics[1])
        networks.append(network)

    return components, networks, auto_deltas

# Construct consensus graph.
def consensus_edges(components, networks):
    edges_to_networks = defaultdict(set)
    for ccs, network in zip(components, networks):
        for cc in ccs:
            for u, v in combinations(cc, 2):
                edges_to_networks[(u, v)].add(network)
    return dict( (edge, len(edge_networks)) for edge, edge_networks in edges_to_networks.iteritems() )
