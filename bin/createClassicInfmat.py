#!/usr/bin/python

import networkx as nx, sys, os, scipy.io
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
    return parser                                                               

def run(args):
    # Load input graph
    print "* Loading input graph..."
    G = nx.Graph()
    G.add_edges_from([map(int, l.rstrip().split()[:2]) for l in open(args.edgelist_file)])
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
    index2gene = hnio.load_index(args.gene_index_file)

    # Compute and save Laplacian
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    print "* Computing Laplacian..."
    L = nx.laplacian_matrix(G)
    scipy.io.savemat("{}/{}_laplacian.mat".format(args.output_dir, args.prefix),
                     dict(L=L),oned_as='column')

    # Exponentiate the Laplacian for the given time and save it
    from scipy.linalg import expm
    Li = expm( -args.time * L )
    scipy.io.savemat("{}/{}_inf_{}.mat".format(args.output_dir, args.prefix, args.time),
                     dict(Li=Li), oned_as='column')

    # Save the index to gene mapping
    index_output_file = "{}/{}_index_genes".format(args.output_dir, args.prefix)
    nodes = G.nodes()
    gene_index_output = ["{} {}".format(i+args.start_index, index2gene[node])
                         for i, node in enumerate(nodes)]
    hnio.write_file(index_output_file, "\n".join(gene_index_output))

    # Create edge list with revised indices
    edge_indices = []
    for u, v in G.edges():
        i = nodes.index(u) + args.start_index
        j = nodes.index(v) + args.start_index
        edge_indices.append( sorted([i, j]) )
    edge_output_file = "{}/{}_edge_list".format(args.output_dir, args.prefix)
    edge_output = ["{} {} 1".format(u, v) for u, v in edge_indices]
    hnio.write_file(edge_output_file, "\n".join(edge_output))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))