#!/usr/bin/python

# Load required modules
import sys, os, numpy as np, networkx as nx, scipy as sp, scipy.io
import os.path
sys.path.append(os.path.split(os.path.split(sys.argv[0])[0])[0])
from hotnet2 import hnap, hnio
from hotnet2.constants import ITERATION_REPLACEMENT_TOKEN

# Parse arguments
def get_parser():
    description = 'Create the personalized PageRank matrix for the given '\
                  'network and restart probability beta.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Location of edgelist file.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Location of gene-index file.')
    parser.add_argument('-pnp', '--permuted_networks_path', required=False, default='',
                        help='Path to influence matrices for permuted networks. Include ' +\
                              ITERATION_REPLACEMENT_TOKEN + ' in the path to be replaced with the\
                              iteration number')
    parser.add_argument('-n', '--network_name', required=True,
                        help='Name of network.')
    parser.add_argument('-o', '--output_file', required=True,
                    help="Output file.")
    parser.add_argument('-s', '--start_index', default=1, type=int,
                    help="Index to output edge list, etc..")
    parser.add_argument('-b', '--beta', required=True, type=float,
                    help="Restart probability beta.")
    parser.add_argument('--exclude_network', required=False, action='store_true',
                    default=False,
                    help="Exclude the gene index and edges in the output file.")
    return parser

# Remove self-loops, multi-edges, and restrict to the largest component
def largest_component(G):
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from( selfLoops )
    return G.subgraph( sorted(nx.connected_components( G ), key=lambda cc: len(cc), reverse=True)[0] )

def run(args):
    # Load gene-index map
    with open(args.gene_index_file) as infile:
        arrs = [ l.rstrip().split() for l in infile ]
        indexToGene = dict((int(arr[0]), arr[1]) for arr in arrs)

    G = nx.Graph()
    G.add_nodes_from( indexToGene.values() ) # in case any nodes have degree zero

    # Load graph
    print "* Loading PPI..."
    with open(args.edgelist_file) as infile:
        edges = [ map(int, l.rstrip().split()[:2]) for l in infile ]
    G.add_edges_from( [(indexToGene[u], indexToGene[v]) for u,v in edges] )

    print "\t- Edges:", len(G.edges())
    print "\t- Nodes:", len(G.nodes())

    # Remove self-loops and restrict to largest connected component
    print "* Removing self-loops, multi-edges, and restricting to",
    print "largest connected component..."
    G = largest_component(G)
    nodes = sorted(G.nodes())
    n = len(nodes)
    print "\t- Largest CC Edges:", len( G.edges() )
    print "\t- Largest CC Nodes:", len( G.nodes() )

    ## Create the PPR matrix either using Scipy or MATLAB
    # Create "walk" matrix (normalized adjacency matrix)
    print "* Creating PPR  matrix..."
    W = nx.to_numpy_matrix( G , nodelist=nodes, dtype=np.float64 )
    W = np.asarray(W)
    W = W / W.sum(axis=0) # normalization step

    ## Create PPR matrix using Python
    from scipy.linalg import inv
    PPR = args.beta*inv(sp.eye(n)-(1.-args.beta)*W)
    if args.exclude_network:
        hnio.save_hdf5(args.output_file, dict(PPR=PPR))
    else:
        output = dict(edges=G.edges(), PPR=PPR, nodes=nodes,
                      network_name=args.network_name,
                      permuted_networks_path=args.permuted_networks_path)
        hnio.save_hdf5(args.output_file, output)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
