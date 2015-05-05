#!/usr/bin/python

import networkx as nx, sys, os, scipy as sp, numpy as np
import os.path
sys.path.append(os.path.split(os.path.split(sys.argv[0])[0])[0])
from hotnet2 import hnap, hnio

def get_parser():                                                             
    description = 'Creates a heat diffusion influence matrix from an input graph.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-o', '--output_dir', required=True, help='Path to output directory.')
    parser.add_argument('-p', '--prefix', required=True, help='Output prefix.')
    parser.add_argument('-s', '--start_index', default=1, type=int,
                        help='Minimum index in the index file.')
    parser.add_argument('-t', '--time', required=True, type=float, help='Diffusion time.')
    parser.add_argument('-f', '--format', default='hdf5', type=str, required=False,
                        choices=['hdf5', 'npy'], help="Output file format.")
    return parser                                                               

def expm_eig(A):
    """
    Compute the matrix exponential for a square, symmetric matrix.
    """
    D, V = sp.linalg.eigh(A)
    return sp.dot(sp.exp(D) * V, sp.linalg.inv(V))

def run(args):
    # Load input graph
    print "* Loading input graph..."
    with open(args.edgelist_file) as infile:
        G = nx.Graph()
        G.add_edges_from([map(int, l.rstrip().split()[:2]) for l in infile])
        print "\t{} nodes with {} edges".format(len(G.nodes()), len(G.edges()))

    # Remove self-loops and zero degree nodes, and
    # restrict to the largest connected component
    print "* Removing self-loops, zero degree nodes, and ",
    print "restricting to the largest connected component"
    G.remove_edges_from([(u,v) for u, v in G.edges() if u == v])
    G.remove_nodes_from([n for n in G.nodes() if G.degree(n) == 0])
    G = G.subgraph(sorted(nx.connected_components( G ), key=lambda cc: len(cc), reverse=True)[0])

    print "\t{} nodes with {} edges remaining".format(len(G.nodes()), len(G.edges()))

    # Load gene index
    indexToGene = hnio.load_index(args.gene_index_file)

    # Compute and save Laplacian
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    print "* Computing Laplacian..."
    L = nx.laplacian_matrix(G)

    # Exponentiate the Laplacian for the given time and save it
    print "* Computing diffusion matrix..."
    Li = expm_eig( -args.time * L.todense() )
    #Li = sp.sparse.linalg.expm( -args.time * L)
    output_prefix = "{}/{}_inf_{}".format(args.output_dir, args.prefix, args.time)
    if args.format == 'hdf5':
        hnio.save_hdf5(output_prefix + ".h5", dict(Li=Li))
    elif args.format == 'npy':
        np.save(output_prefix + ".npy", Li)

    # Save the index to gene mapping
    indexOutputFile = "{}/{}_index_genes".format(args.output_dir, args.prefix)
    nodes = G.nodes()
    geneIndexOutput = ["{} {}".format(i+args.start_index, indexToGene[node])
                         for i, node in enumerate(nodes)]
    hnio.write_file(indexOutputFile, "\n".join(geneIndexOutput))

    # Create edge list with revised indices
    edgeIndices = []
    for u, v in G.edges():
        i = nodes.index(u) + args.start_index
        j = nodes.index(v) + args.start_index
        edgeIndices.append( sorted([i, j]) )
    edgeOutputFile = "{}/{}_edge_list".format(args.output_dir, args.prefix)
    edgeOutput = ["{} {} 1".format(u, v) for u, v in edgeIndices]
    hnio.write_file(edgeOutputFile, "\n".join(edgeOutput))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
