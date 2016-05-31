#!/usr/bin/python

# Load required modules
import sys, os, numpy as np, networkx as nx, scipy as sp, scipy.io
import os.path
sys.path.append(os.path.split(os.path.split(sys.argv[0])[0])[0])
from hotnet2 import hnap, hnio

# Parse arguments
def get_parser():
    description = 'Create the personalized PageRank matrix for the given '\
                  'network and restart probability beta.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Location of edgelist file.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Location of gene-index file.')
    parser.add_argument('-o', '--output_dir', required=True,
	                help="Output dir.")
    parser.add_argument('-p', '--prefix', required=True,
	                help="Output prefix.")
    parser.add_argument('-s', '--start_index', default=1, type=int,
	                help="Index to output edge list, etc..")
    parser.add_argument('-b', '--beta', required=True, type=float,
	                help="Restart probability beta.")
    parser.add_argument('-f', '--format', default='hdf5', type=str,
                    choices=['hdf5', 'npy', 'matlab'],
                    help="Output file format.")
    return parser

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
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from( selfLoops )
    G = G.subgraph( sorted(nx.connected_components( G ), key=lambda cc: len(cc),
                           reverse=True)[0] )
    nodes = sorted(G.nodes())
    n = len(nodes)
    print "\t- Largest CC Edges:", len( G.edges() )
    print "\t- Largest CC Nodes:", len( G.nodes() )

    # Set up output directory
    print "* Saving updated graph to file..."
    os.system( 'mkdir -p ' + args.output_dir )
    output_dir = os.path.normpath(args.output_dir)
    output_prefix = "{}/{}".format(output_dir, args.prefix)

    if args.format == 'hdf5': ext = 'h5'
    elif args.format == 'matlab': ext = 'mat'
    else: ext = args.format

    pprfile = "{}_ppr_{:g}.{}".format(output_prefix, args.beta, ext)

    # Index mapping for genes
    index_map = [ "{} {}".format(i+args.start_index, nodes[i]) for i in range(n) ]
    with open("{}_index_genes".format(output_prefix), 'w') as outfile:
        outfile.write( "\n".join(index_map) )

    # Edge list
    edges = [sorted([nodes.index(u) + args.start_index,
                     nodes.index(v) + args.start_index])
             for u, v in G.edges()]
    edgelist = [ "{} {} 1".format(u, v) for u, v in edges ]

    with open("{}_edge_list".format(output_prefix), 'w') as outfile:
        outfile.write( "\n".join(edgelist) )

    ## Create the PPR matrix either using Scipy or MATLAB
    # Create "walk" matrix (normalized adjacency matrix)
    print "* Creating PPR  matrix..."
    W = nx.to_numpy_matrix( G , nodelist=nodes, dtype=np.float64 )
    W = np.asarray(W)
    W = W / W.sum(axis=0) # normalization step

    ## Create PPR matrix using Python
    from scipy.linalg import inv
    PPR = args.beta*inv(sp.eye(n)-(1.-args.beta)*W)
    if args.format == 'hdf5':
        hnio.save_hdf5(pprfile, dict(PPR=PPR))
    elif args.format == 'npy':
        np.save(pprfile, PPR)
    elif args.format == 'matlab':
        scipy.io.savemat(pprfile, dict(PPR=PPR))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
