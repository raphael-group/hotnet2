#!/usr/bin/env python

import numpy as np, random, networkx as nx
from itertools import combinations
import os, sys, argparse

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_vertices', type=int, required=True)
    parser.add_argument('-m', '--min_num_edges', type=int, required=True)
    parser.add_argument('-i', '--implant_sizes', type=int, required=True, nargs='*')
    parser.add_argument('-s', '--seed', type=int, required=False, default=0)
    parser.add_argument('-nsn', '--num_sampled_vertices_network', type=int, nargs='*', required=True)
    parser.add_argument('-nsw', '--num_sampled_vertices_weights', type=int, nargs='*', required=True)
    parser.add_argument('-ivf', '--index_vertex_files', type=str, nargs='*', required=True)
    parser.add_argument('-elf', '--edge_list_files', type=str, nargs='*', required=True)
    parser.add_argument('-vwf', '--vertex_weight_files', type=str, nargs='*', required=True)
    return parser

def letter_range(*args):
    """
    Return range with (combinations of) letters instead of integers.

    Arguments (identical for range; depends on number of given arguments):
    m: start
    n: stop
    s: step

    Example usage:
    In [1]: letter_range(5)
    Out[1]: ['a', 'b', 'c', 'd', 'e']
    In [2]: letter_range(24, 28)
    Out[2]: ['y', 'z', 'aa', 'ab']
    In [3]: letter_range(0, 28, 5)
    Out[3]: ['a', 'f', 'k', 'p', 'u', 'z']
    In [4]: letter_range(12345)[-6:]
    Out[4]: ['cssx', 'cssy', 'cssz', 'cstt', 'cstu', 'cstv']
    """
    import itertools, string

    if len(args)==1:
        m, n, s = 0, args[0], 1
    elif len(args)==2:
        m, n, s = args[0], args[1], 1
    elif len(args)==3:
        m, n, s = args[0], args[1], args[2]
    else:
        raise TypeError('letter_range expects 1, 2, or 3 arguments but received {} arguments.'.format(len(args)))

    letters = list(string.ascii_lowercase)
    result = []

    i = 0
    for k in range(1, n+1):
        for subset in itertools.combinations_with_replacement(letters, k):
            if i>=m and (s==1 or i%s==0):
                result.append(''.join(subset))
            i += 1
            if i==n:
                return result

def run(args):
    # Check arguments.
    assert len(args.num_sampled_vertices_network)==len(args.index_vertex_files)==len(args.edge_list_files)
    assert len(args.num_sampled_vertices_weights)==len(args.vertex_weight_files)
    assert args.num_vertices>0
    assert args.min_num_edges>0
    assert args.min_num_edges<=args.num_vertices*(args.num_vertices-1)/2

    # Create random network.
    n = args.num_vertices
    for i in range(1, n):
        G = nx.barabasi_albert_graph(n, i, seed=args.seed)
        if nx.number_of_edges(G)>args.min_num_edges:
            break

    # Use letters instead of integers for the vertex labels.
    vertices = letter_range(n)
    edges = [(vertices[i], vertices[j]) for i, j in G.edges()]
    G = nx.Graph()
    G.add_nodes_from(vertices)
    G.add_edges_from(edges)

    # Implant cliques.
    random.seed(args.seed)
    implanted_vertices = set()
    nonimplanted_vertices = set(vertices)
    for implant_size in args.implant_sizes:
        implant = set(random.sample(sorted(nonimplanted_vertices), implant_size))
        implanted_vertices = set.union(implanted_vertices, implant)
        nonimplanted_vertices = set.difference(nonimplanted_vertices, implant)
        for u, v in combinations(implant, 2):
            G.add_edge(u, v)
    implanted_vertices = sorted(implanted_vertices)
    nonimplanted_vertices = sorted(nonimplanted_vertices)

    # Create random weights.
    np.random.seed(args.seed)
    weights = np.random.exponential(scale=1.0, size=n)

    # Assign high weights to implanted vertices and low weights to other vertices.
    random.seed(args.seed)
    sorted_weights = np.sort(weights)[::-1]
    random.shuffle(implanted_vertices)
    random.shuffle(nonimplanted_vertices)
    sorted_vertices = implanted_vertices + nonimplanted_vertices
    vertex_to_weight = dict((v, w) for v, w in zip(sorted_vertices, sorted_weights))

    # Sample the random network.
    G_samples = []
    for num_sampled_vertices in args.num_sampled_vertices_network:
        random.seed(args.seed)
        sampled_vertices = random.sample(vertices, num_sampled_vertices)
        G_sample = G.subgraph(sampled_vertices)
        G_samples.append(G_sample)

    # Sample the random weights.
    vertex_to_weight_samples = []
    for num_sampled_vertices in args.num_sampled_vertices_weights:
        random.seed(args.seed)
        sampled_vertices = random.sample(vertices, num_sampled_vertices)
        vertex_to_weight_sample = dict((v, vertex_to_weight[v]) for v in sampled_vertices)
        vertex_to_weight_samples.append(vertex_to_weight_sample)

    # Output sampled random networks.
    for G_sample, index_vertex_file, edge_list_file in zip(G_samples, args.index_vertex_files, args.edge_list_files):
        np.random.seed(args.seed)
        sampled_vertices = np.random.permutation(G_sample.nodes())
        vertex_to_index = dict((v, i+1) for i, v in enumerate(sampled_vertices))
        edge_list = [[vertex_to_index[u], vertex_to_index[v]] for u, v in G_sample.edges() if u!=v]
        with open(index_vertex_file, 'w') as f:
            index_vertex_string = '\n'.join('\t'.join(map(str, (i, v))) for v, i in vertex_to_index.iteritems())
            f.write(index_vertex_string)
        with open(edge_list_file, 'w') as f:
            edge_list_string = '\n'.join('\t'.join(map(str, edge)) for edge in edge_list)
            f.write(edge_list_string)

    # Output sample random weights.
    for vertex_to_weight_sample, vertex_weight_file in zip(vertex_to_weight_samples, args.vertex_weight_files):
        with open(vertex_weight_file, 'w') as f:
            vertex_weight_string = '\n'.join('\t'.join(map(str, (v, w))) for v, w in vertex_to_weight_sample.iteritems())
            f.write(vertex_weight_string)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
