# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy as sp
import scipy.io
import hotnet2 as hn
import networkx
strong_ccs = networkx.strongly_connected_components

def create_permuted_sim_mat(permuted_mat, permuted_mat_index, genes, h):
    M, gene_index = hn.induce_infmat(permuted_mat, permuted_mat_index, genes)
    sim_mat = hn.similarity_matrix(M, h, gene_index)

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
        delta = -1
        index = round(len(sorted_edges)* start_quant)
        left, right = 0., float(len(sorted_edges))      #TODO: why are these floats rather than ints?
        visited = []
        
        while len(visited) < 100:
            # print "\t\t\t%s: %s" % (len(visited), index/len(sorted_edges)),
            # print "(%s)" % format(sorted_edges[int(index)], 'g')
            
            # construct graph using new delta
            delta = sorted_edges[int(index)]
            DiG = hn.weighted_graph(permuted_sim, permuted_index, delta)
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
    delta_maps = map_fn(network_delta_wrapper, args)
    
    if parallel:
        pool.close()
        pool.join()
        
    # Parse the delta_maps into one dictionary
    sizes2deltas = dict([(s, []) for s in sizes])
    for size2delta in delta_maps:
        for s in sizes: sizes2deltas[s].append( size2delta[s] )

    return sizes2deltas


def heat_delta_wrapper_old( (M, h, gene_index, sizes, component_fn) ):
    permuted_h = np.array([ val for val in h] )
    shuffle( permuted_h )
    sim_mat = hn.similarity_matrix( M, permuted_h, gene_index )
    return find_best_delta( sim_mat, gene_index, sizes, component_fn )

def heat_delta_wrapper((M, heat_permutation)):
    pass

def heat_delta_selection(M, gene_index, heat_permutations, sizes,
                         component_fn=strong_ccs, parallel=True):
    h = hn.heat_vec(gene2heat, gene_index)
    return None

#list of num_permutations dicts of max cc size => best delta
from random import shuffle
def heat_delta_selection_old(M, gene_index, h, num_permutations, sizes,
                             component_fn=strong_ccs, parallel=True):
    print "* Performing permuted heat delta selection..."
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map

    args = [(M, h, gene_index, sizes, component_fn)] * num_permutations
    deltas = map_fn(heat_delta_wrapper_old, args)

    if parallel:
        pool.close()
        pool.join()
        
    return deltas
