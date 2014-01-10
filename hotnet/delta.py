# -*- coding: iso-8859-1 -*-
import scipy as sp
import scipy.io
import numpy as np
import hotnet2 as hn
import networkx as nx
import multiprocessing as mp
strong_ccs = nx.strongly_connected_components
from union_find import UnionFind
from collections import namedtuple, defaultdict

def get_component_sizes(arrs):
    return [len(arr) for arr in arrs]

def delta_too_small(component_sizes, max_size ):
    # print "\t\t\t", max(component_sizes)
    if max_size > max(component_sizes): return True
    else: return False

def find_best_delta_by_largest_cc(permuted_sim, permuted_index, sizes, start_quant=0.99):
    """Return a dict mapping each size in sizes to the smallest delta such that the size of the
    largest CC in the graph corresponding the the given similarity matrix is <= that size.
    
    Arguments:
    permuted_sim -- 2D ndarray representing a similarity matrix
    permuted_index -- dict mapping an index in the matrix to the name of the gene represented at that
                      index in the influence matrix
    sizes -- list of sizes for the largest connected component
    start_quant -- percentile of edge weights that should be used as the starting delta for the
                   binary search procedure
    
    """
    
    print "Finding smallest delta such that size of largest CC is <= l"
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
            G = hn.weighted_graph(permuted_sim, permuted_index, delta)
            if delta in visited:
                size2delta[max_size] = delta
                break
            else:
                visited.append( delta )
            
            # increment / decrement index based on whether components meet size specifications
            if delta_too_small( get_component_sizes( strong_ccs( G) ), max_size ):
                right  = index
                index -= round( (index-left)/2. )
            else:
                left   = index
                index += round( (right-index)/2. )

        if max_size not in size2delta:
            raise AssertionError("NO DELTA SELECTED FOR k = " + str(max_size))
            
    return size2delta

Edge = namedtuple("Edge", ["node1", "node2", "weight"])

def network_delta_wrapper((network_path, infmat_name, index2gene, heat, sizes, selection_function)):
    permuted_mat = scipy.io.loadmat(network_path)[infmat_name]   
    M, gene_index = hn.induce_infmat(permuted_mat, index2gene, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h)
    if selection_function is find_best_delta_by_largest_cc:
        return selection_function(sim, gene_index, sizes)
    else:
        raise ValueError("Unknown delta selection function: %s" % (selection_function))

def network_delta_selection(network_paths, infmat_name, index2gene, heat, sizes,
                            parallel=True, selection_fn=find_best_delta_by_largest_cc):
    """Return a dict mapping each size in sizes to a list of the best deltas for each permuted
    network for that size. 
    
    Arguments:
    network_paths -- iterable of paths to .mat files containing permuted networks
    infmat_name -- name of influence matrix in .mat files
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    heat -- dict mapping a gene name to the heat score for that gene
    sizes -- list of sizes for largest CC / min size for CCs to be counted (based on selection_fn)
    parallel -- whether finding the best delta for each permuted network should be performed in parallel
    selection_fn -- function that should be used for finding the best delta 
    
    """
    if parallel:
        pool = mp.Pool(25)
        map_fn = pool.map
    else:
        map_fn = map
       
    args = [(network_path, infmat_name, index2gene, heat, sizes, selection_fn)
            for network_path in network_paths]
    delta_maps = map_fn(network_delta_wrapper, args)
    
    if parallel:
        pool.close()
        pool.join()
         
    # Parse the deltas into one dictionary
    sizes2deltas = dict([(s, []) for s in sizes])
    for size2delta in delta_maps:
        for s in sizes: sizes2deltas[s].append( size2delta[s] )
 
    return sizes2deltas
