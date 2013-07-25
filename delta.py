# -*- coding: iso-8859-1 -*-
import scipy as sp
import scipy.io
import numpy as np
import hotnet2 as hn
import networkx as nx
import multiprocessing as mp
strong_ccs = nx.strongly_connected_components
from union_find import UnionFind

def get_component_sizes(arrs):
    return [len(arr) for arr in arrs]

def delta_too_small(component_sizes, max_size ):
    # print "\t\t\t", max(component_sizes)
    if max_size > max(component_sizes): return True
    else: return False

def find_best_delta_by_largest_cc(permuted_sim, permuted_index, sizes, directed, start_quant=0.99):
    component_fn = strong_ccs if directed else nx.connected_components
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

def find_best_delta_by_num_ccs(permuted_sim, k, start=0.05):
    if k < 2:
        raise ValueError("k must be at least 2")
    
    edges = get_edges(permuted_sim, start)
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
    """Return a list of Edge objects representing edges in the top start% of edge weights"""
    flattened = np.ndarray.flatten(sim)
    if np.array_equal(sim, sim.transpose()):
        edges = [Edge(i/len(sim), i%len(sim), flattened[i]) for i in range(len(flattened)) if i/len(sim) <= i%len(sim)]
    else:
        edges = [Edge(i/len(sim), i%len(sim), flattened[i]) for i in range(len(flattened))]
    edges = sorted(edges, key=lambda x: x.weight, reverse=True)
    edges = edges[:int(start*len(edges))]
    return edges

class Edge:
    def __init__(self, node1, node2, weight):
        self.node1 = node1
        self.node2 = node2
        self.weight = weight

    def __repr__(self):
        return '(%s, %s, %s)' % (self.node1, self.node2, self.weight)

def network_delta_wrapper((network_path, infmat_name, index2gene, tested_genes, h, sizes, directed)):
    permuted_mat = scipy.io.loadmat(network_path)[infmat_name]   
    M, gene_index = hn.induce_infmat(permuted_mat, index2gene, tested_genes)
    sim = hn.similarity_matrix(M, h, directed)
    return find_best_delta_by_largest_cc(sim, gene_index, sizes, directed )

def network_delta_selection(network_paths, infmat_name, index2gene, heat, sizes,
                            directed=True, parallel=True):
    print "* Performing network delta selection..."
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map
       
    h_vec = hn.heat_vec(heat, index2gene)
    args = [(network_path, infmat_name, index2gene, sorted(heat.keys()), h_vec, sizes, directed)
            for network_path in network_paths]
    delta_maps = map_fn(network_delta_wrapper, args)
    
    if parallel:
        pool.close()
        pool.join()
         
    # Parse the delta_maps into one dictionary
    sizes2deltas = dict([(s, []) for s in sizes])
    for size2delta in delta_maps:
        for s in sizes: sizes2deltas[s].append( size2delta[s] )
 
    return sizes2deltas


def heat_delta_wrapper((M, index2gene, heat_permutation, directed, sizes, selection_function)):
    heat = hn.heat_vec(heat_permutation, index2gene)
    sim_mat = hn.similarity_matrix(M, heat, directed)
    return selection_function(sim_mat, index2gene, sizes, directed)

#list of num_permutations dicts of max cc size => best delta
def heat_delta_selection(M, index2gene, heat_permutations, sizes, directed=True, parallel=True,
                         selection_fn=find_best_delta_by_largest_cc):
    print "* Performing permuted heat delta selection..."
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map

    args = [(M, index2gene, heat_permutation, directed, sizes, selection_fn)
            for heat_permutation in heat_permutations]
    deltas = map_fn(heat_delta_wrapper, args)

    if parallel:
        pool.close()
        pool.join()
        
    return deltas
