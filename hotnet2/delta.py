# -*- coding: iso-8859-1 -*-
from collections import namedtuple, defaultdict
import multiprocessing as mp
import numpy as np
import scipy as sp
import networkx as nx
import hotnet2 as hn
import hnio
from union_find import UnionFind

strong_ccs = nx.strongly_connected_components

def get_component_sizes(arrs):
    return [len(arr) for arr in arrs]

def delta_too_small(component_sizes, max_size ):
    # print "\t\t\t", max(component_sizes)
    if max_size > max(component_sizes): return True
    else: return False

def find_best_delta_by_largest_cc(permuted_sim, permuted_index, sizes, directed, start_quant=0.99):
    """Return a dict mapping each size in sizes to the smallest delta such that the size of the
    largest CC in the graph corresponding the the given similarity matrix is <= that size.
    
    Arguments:
    permuted_sim -- 2D ndarray representing a similarity matrix
    permuted_index -- dict mapping an index in the matrix to the name of the gene represented at that
                      index in the influence matrix
    sizes -- list of sizes for the largest connected component
    directed -- whether or not the graph constructed from the similarity matrix should be directed
    start_quant -- percentile of edge weights that should be used as the starting delta for the
                   binary search procedure
    
    """
    
    print "Finding smallest delta such that size of largest CC is <= l"
    component_fn = strong_ccs if directed else nx.connected_components
    # Construct weighted digraphs for each network for each delta
    sorted_edges = np.unique(permuted_sim) # unique edge weights in sorted ascending
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
            G = hn.weighted_graph(permuted_sim, permuted_index, delta, directed)
            if delta in visited:
                size2delta[max_size] = delta
                break
            else:
                visited.append( delta )
            
            # increment / decrement index based on whether components meet size specifications
            if delta_too_small( get_component_sizes( component_fn( G) ), max_size ):
                right  = index
                index -= round( (index-left)/2. )
            else:
                left   = index
                index += round( (right-index)/2. )

        if max_size not in size2delta:
            raise AssertionError("NO DELTA SELECTED FOR k = " + str(max_size))
            
    return size2delta

Edge = namedtuple("Edge", ["node1", "node2", "weight"])

def find_best_delta_by_num_ccs(permuted_sim, ks, start=0.05):
    """Return a dict mapping each size in ks to the median delta that maximizes the number of
    connected components of size at least k in the graph corresponding the the given similarity
    matrix.
    
    Arguments:
    permuted_sim -- 2D ndarray representing a similarity matrix
    ks -- list of minimum sizes for connected components to be counted. This must be at least 2.
    start -- only deltas in the top start proportion of edge weights will be considered when
             searching for deltas that maximize the number of connected components
    """
    
    print "Finding median delta that maximizes the # of CCs of size >= l"
    edges = get_edges(permuted_sim, start)
    k2delta = {}

    for k in ks:
        _, bestDeltas = find_best_delta_by_num_ccs_for_given_k(permuted_sim, edges, k)
        k2delta[k] = np.median(bestDeltas)

    return k2delta

def find_best_delta_by_num_ccs_for_given_k(permuted_sim, edges, k):
    
    if k < 2:
            raise ValueError("k must be at least 2")

    max_num_ccs = 0 #initially, each node is its own CC of size 1, so none is of size >= k for k >= 2
    bestDeltas = [edges[0].weight]
    uf = UnionFind()

    for edge in edges:
        uf.union(edge.node1, edge.node2)
        num_ccs = len([root for root in uf.roots if uf.weights[root] >= k])
        if num_ccs > max_num_ccs:
            max_num_ccs = num_ccs
            bestDeltas = [edge.weight]
        elif num_ccs == max_num_ccs:
            bestDeltas.append(edge.weight)

    return max_num_ccs, bestDeltas

def get_edges(sim, start=.05):
    """Return a list of Edge tuples representing edges in the top start% of edge weights"""
    flattened = np.ndarray.flatten(sim)
    if np.array_equal(sim, sim.transpose()):
        edges = [Edge(i/len(sim), i%len(sim), flattened[i]) for i in range(len(flattened)) if i/len(sim) <= i%len(sim)]
    else:
        edges = [Edge(i/len(sim), i%len(sim), flattened[i]) for i in range(len(flattened))]
    edges = sorted(edges, key=lambda x: x.weight, reverse=True)
    edges = edges[:int(start*len(edges))]
    return edges

def network_delta_wrapper((network_path, infmat_name, index2gene, heat, sizes, directed,
                           selection_function)):
    permuted_mat = hnio.load_infmat(network_path, infmat_name)
    sim, index2gene = hn.similarity_matrix(permuted_mat, index2gene, heat, directed)
    if selection_function is find_best_delta_by_largest_cc:
        return selection_function(sim, index2gene, sizes, directed)
    elif selection_function is find_best_delta_by_num_ccs:
        return selection_function(sim, sizes)
    else:
        raise ValueError("Unknown delta selection function: %s" % (selection_function))

def network_delta_selection(network_paths, infmat_name, index2gene, heat, sizes, directed=True,
                            num_cores=1, selection_fn=find_best_delta_by_largest_cc):
    """Return a dict mapping each size in sizes to a list of the best deltas for each permuted
    network for that size. 
    
    Arguments:
    network_paths -- iterable of paths to .mat files containing permuted networks
    infmat_name -- name of influence matrix in .mat files
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    heat -- dict mapping a gene name to the heat score for that gene
    sizes -- list of sizes for largest CC / min size for CCs to be counted (based on selection_fn)
    directed -- whether or not the graph constructed from the similarity matrix should be directed
    num_cores -- number of cores to use for running in parallel
    selection_fn -- function that should be used for finding the best delta 
    
    """
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map
       
    args = [(network_path, infmat_name, index2gene, heat, sizes, directed,
             selection_fn) for network_path in network_paths]
    delta_maps = map_fn(network_delta_wrapper, args)
    
    if num_cores != 1:
        pool.close()
        pool.join()
         
    # Parse the deltas into one dictionary
    sizes2deltas = defaultdict(list)
    for size2delta in delta_maps:
        for s in sizes: sizes2deltas[s].append(size2delta[s])
 
    return sizes2deltas

def heat_delta_wrapper((infmat, index2gene, heat_permutation, directed, sizes, selection_function)):
    sim, index2gene = hn.similarity_matrix(infmat, index2gene, heat_permutation, directed)
    if selection_function is find_best_delta_by_largest_cc:
        return selection_function(sim, index2gene, sizes, directed)
    elif selection_function is find_best_delta_by_num_ccs:
        return selection_function(sim, sizes)
    else:
        raise ValueError("Unknown delta selection function: %s" % (selection_function))

#list of num_permutations dicts of max cc size => best delta
def heat_delta_selection(infmat, index2gene, heat_permutations, sizes, directed=True, num_cores=1,
                         selection_fn=find_best_delta_by_largest_cc):
    """Return a dict mapping each size in sizes to a list of the best deltas for each heat
    permutation for that size. 
    
    Arguments:
    infmat -- 2D ndarray representing an influence matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    heat_permutations -- list of heat permutations (dicts mapping gene name to heat score)
    sizes -- list of sizes for largest CC / min size for CCs to be counted (based on selection_fn)
    directed -- whether or not the graph constructed from the similarity matrix should be directed
    num_cores -- number of cores to use for running in parallel
    selection_fn -- function that should be used for finding the best delta 
    
    """
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map
    args = [(infmat, index2gene, heat_permutation, directed, sizes, selection_fn)
            for heat_permutation in heat_permutations]
    deltas = map_fn(heat_delta_wrapper, args)

    if num_cores != 1:
        pool.close()
        pool.join()
    
    # Parse the deltas into one dictionary
    sizes2deltas = defaultdict(list)
    for size2delta in deltas:
        for s in sizes: sizes2deltas[s].append(size2delta[s])
         
    return sizes2deltas
