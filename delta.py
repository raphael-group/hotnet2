#!/usr/bin/python

import numpy as np, sys
from data import load_index
from hotnet2 import *

## Data locations ##
# Interaction network directories
KERNEL_DIR   = "/data/compbio/datasets/HeatKernels/pagerank/permuted_networks/pprs"
IREF_DIR     = KERNEL_DIR + "/iref/"
HINT_DIR     = KERNEL_DIR + "/inthint/"
MULTINET_DIR = KERNEL_DIR + "/multinet/"

# influence matrix locations
pprmats       = { IREF: {}, HINT: {}, MULTINET: {} }
ppr_edgelists = { IREF: {}, HINT: {}, MULTINET: {} }
ppr_indices   = { IREF: {}, HINT: {}, MULTINET: {} }
alphas = ['%.2f' % alpha for alpha in np.arange(0.05, 0.95, 0.05)]
for alpha in alphas:
    for network in [IREF, HINT, MULTINET]:
        pprmats[network][alpha], ppr_edgelists[network][alpha] = {}, {}
        ppr_indices[network][alpha] = {}
        
for alpha in alphas:
    for i in range(1, 101):
        pprmats[IREF][alpha][i] = "%s/%s/iref_ppr_%s.mat" % (IREF_DIR, i, alpha)
        ppr_edgelists[IREF][alpha][i] = "%s/%s/iref_edge_list" % (IREF_DIR, i)
        ppr_indices[IREF][alpha][i] = "%s/%s/iref_index_genes" % (IREF_DIR, i)
        pprmats[HINT][alpha][i] = "%s/%s/inthint_ppr_%s.mat" % (HINT_DIR, i, alpha)
        ppr_edgelists[HINT][alpha][i] = "%s/%s/inthint_edge_list" % (HINT_DIR, i)
        ppr_indices[HINT][alpha][i] = "%s/%s/inthint_index_genes" % (HINT_DIR, i)
        pprmats[MULTINET][alpha][i] = "%s/%s/multinet_ppr_%s.mat" % (MULTINET_DIR, i, alpha)
        ppr_edgelists[MULTINET][alpha][i] = "%s/%s/multinet_edge_list" % (MULTINET_DIR, i)
        ppr_indices[MULTINET][alpha][i] = "%s/%s/multinet_index_genes" % (MULTINET_DIR, i)


## Manipulate permuted infmats ##
def load_permuted_pprmat(network, beta, num):
    alpha = '%.2f' % (1-beta)
    infmat = scipy.io.loadmat( pprmats[network][alpha][num] )["PPR"]
    gene_index = load_index( ppr_indices[network][alpha][num] )
    return gene_index, infmat

def create_permuted_sim_mat( permuted_mat, permuted_mat_index, h, genes ):
    # Create an induced similarity matrix
    M, gene_index, _ = induce_infmat( permuted_mat, permuted_mat_index, genes )
    sim_mat, _ = similarity_matrix( M, h, gene_index )

    return sim_mat, gene_index

def get_component_sizes(arrs): return [len(arr) for arr in arrs]

def delta_too_small(component_sizes, max_size ):
    # print "\t\t\t", max(component_sizes)
    if max_size > max(component_sizes): return True
    else: return False

def find_best_delta(permuted_sim, permuted_index, sizes, component_fn, start_quant=0.99):
    # Construct weighted digraphs for each network for each delta
    sorted_edges = sorted(sp.ndarray.flatten(permuted_sim))
    size2delta = dict()
    for max_size in sizes:
        # print "\t\tk=", max_size
        prev_delta, delta = -1, -1
        index = round(len(sorted_edges)* start_quant)
        left, right = 0., float(len(sorted_edges))
        visited = []
        
        while len(visited) < 100:
            # print "\t\t\t%s: %s" % (len(visited), index/len(sorted_edges)),
            # print "(%s)" % format(sorted_edges[int(index)], 'g')
            prev_index = index
            
            # construct graph using new delta
            delta = sorted_edges[int(index)]
            DiG   = weighted_graph( permuted_sim, permuted_index, delta )
            if delta in visited:
                size2delta[max_size] = delta
                break
            else:
                visited.append( delta )
            
            # increment / decrement index based on whether components meet size specifications
            if delta_too_small( get_component_sizes( component_fn( DiG) ), max_size ):
                right  = index
                index -= round( (index-left)/2. )
            else:
                left   = index
                index += round( (right-index)/2. )

        if max_size not in size2delta:
            print "NO DELTA SELECTED FOR k =", max_size
            sys.exit()
            
    return size2delta

def network_delta_wrapper( (network, beta, h, genes, i, sizes, component_fn) ):
    permuted_mat_index, permuted_mat = load_permuted_pprmat( network, beta, i )
    sim_mat, gene_index = create_permuted_sim_mat(permuted_mat, permuted_mat_index, h, genes)
    return find_best_delta( sim_mat, gene_index, sizes, component_fn )

import multiprocessing as mp
def network_delta_selection( network, beta, h, genes, num_permutations, sizes,
                             component_fn=strong_ccs, parallel=True):
    print "* Performing network delta selection..."
    if parallel:
        # Find the optimal delta for each of the given input sizes
        pool = mp.Pool()
        args = [(network, beta, h, genes, i+1, sizes, component_fn)
                for i in range(num_permutations)]
        delta_maps = pool.map( network_delta_wrapper, args )
        pool.close()
        pool.join()

        # Parse the delta_maps into one dictionary
        sizes2deltas = dict([(s, []) for s in sizes])
        for size2delta in delta_maps:
            for s in sizes: sizes2deltas[s].append( size2delta[s] )

        return sizes2deltas

    else:
        sizes2deltas = dict([(size, []) for size in sizes])
        for i in range(num_permutations):
            print "\t- Permutation", i+1
            permuted_mat_index, permuted_mat = load_permuted_pprmat( network, beta, i+1 )
            sim_mat, gene_index = create_permuted_sim_mat(permuted_mat, permuted_mat_index, h,
                                                          genes)
            size2delta = find_best_delta( sim_mat, gene_index, sizes, component_fn)
            print "\t\t=>Best deltas:"
            for size in sizes:
                print "\t\t\tk=%s: %s" % (size, size2delta[size])
                sizes2deltas[size].append( size2delta[size] )
        
    return sizes2deltas


def heat_delta_wrapper( (M, h, gene_index, component_fn, sizes) ):
    permuted_h = np.array([ val for val in h] )
    shuffle( permuted_h )
    sim_mat, _ = similarity_matrix( M, permuted_h, gene_index )
    return find_best_delta( sim_mat, gene_index, sizes, component_fn )

from random import shuffle
def heat_delta_selection( M, gene_index, h, num_permutations, sizes,
                          component_fn=strong_ccs, parallel=True):
    print "* Performing permuted heat delta selection..."
    if parallel:
        pool = mp.Pool()
        args = [(M, h, gene_index, sizes, component_fn )
                for i in range(num_permutations)]
        deltas = pool.map( heat_delta_wrapper, args )
        pool.close()
        pool.join()

    else:
        deltas = []
        for i in range(1, num_permutations+1):
            print "\t- Permutation", i
            permuted_h = np.array([ val for val in h] )
            shuffle( permuted_h )
            sim_mat, _ = similarity_matrix( M, permuted_h, gene_index )
            best_delta = find_best_delta( sim_mat, gene_index, sizes, component_fn )
            deltas.append( best_delta )
            print "\t\t=>Best delta:", best_delta
        
    return deltas
