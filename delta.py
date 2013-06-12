# -*- coding: iso-8859-1 -*-

import numpy as np, sys
import scipy.io
from hotnet2 import *

def create_permuted_sim_mat(permuted_mat, permuted_mat_index, genes, h):
    M, gene_index, _ = induce_infmat(permuted_mat, permuted_mat_index, genes)
    sim_mat, _ = similarity_matrix(M, h, gene_index)

    return sim_mat, gene_index

def get_component_sizes(arrs):
    return [len(arr) for arr in arrs]

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
            raise AssertionError("NO DELTA SELECTED FOR k = " + str(max_size))
            
    return size2delta


def network_delta_wrapper((network_path, infmat_name, index2gene, tested_genes, h, sizes, component_fn)):
    permuted_mat = scipy.io.loadmat(network_path)[infmat_name]    
    sim_mat, gene_index = create_permuted_sim_mat(permuted_mat, index2gene, tested_genes, h)
    return find_best_delta( sim_mat, gene_index, sizes, component_fn )


import multiprocessing as mp
def network_delta_selection(permuted_network_paths, index2gene, infmat_name, tested_genes, h, sizes,
                            component_fn=strong_ccs, parallel=True):
    print "* Performing network delta selection..."
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map
        
    args = [(network_path, infmat_name, index2gene, tested_genes, h, sizes, component_fn) for network_path in permuted_network_paths]
    delta_maps = map_fn( network_delta_wrapper, args )
    
    if parallel:
        pool.close()
        pool.join()
        
    # Parse the delta_maps into one dictionary
    sizes2deltas = dict([(s, []) for s in sizes])
    for size2delta in delta_maps:
        for s in sizes: sizes2deltas[s].append( size2delta[s] )

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
