#!/usr/bin/env python

# Load required modules
import networkx as nx, numpy as np, h5py
from scipy.linalg import inv
import hnio

# Run the entire HotNet2 diffusion process from start to finish
def save_hotnet2_diffusion_to_file( index_file, edge_file, beta, output_file,
                                    exclude_network=False, params=dict(), verbose=0):
    # Load the graph
    if verbose > 0: print "* Loading PPI..."
    G = load_network_from_file( index_file, edge_file)
    if verbose > 0:
        print "\t- Edges:", len(G.edges())
        print "\t- Nodes:", len(G.nodes())

    # Remove self-loops and restrict to largest connected component
    if verbose > 0:
        print "* Removing self-loops, multi-edges, and restricting to",
        print "largest connected component..."

    G = largest_component(G)
    nodes = sorted(G.nodes())
    n = len(nodes)

    if verbose > 0:
        print "\t- Largest CC Edges:", len( G.edges() )
        print "\t- Largest CC Nodes:", len( G.nodes() )

    # Run the diffusion and output to file
    PPR = hotnet2_diffusion(G, nodes, beta, verbose)

    if exclude_network:
        hnio.save_hdf5(output_file, dict(PPR=PPR))
    else:
        output = dict(edges=G.edges(), PPR=PPR, nodes=nodes)
        output.update(params.items())
        hnio.save_hdf5(output_file, output)

# Load a graph from file
def load_network_from_file(index_file, edge_file):
    # Load gene-index map
    with open(index_file) as infile:
        arrs = [ l.rstrip().split() for l in infile ]
        indexToGene = dict((int(arr[0]), arr[1]) for arr in arrs)

    G = nx.Graph()
    G.add_nodes_from( indexToGene.values() ) # in case any nodes have degree zero

    # Load graph
    with open(edge_file) as infile:
        edges = [ map(int, l.rstrip().split()[:2]) for l in infile ]
    G.add_edges_from( [(indexToGene[u], indexToGene[v]) for u,v in edges] )

    return G

# Run the HotNet2 diffusion process on a given network
def hotnet2_diffusion(G, nodes, beta, verbose):
    if verbose > 1: print "* Creating PPR  matrix..."
    W = nx.to_numpy_matrix( G , nodelist=nodes, dtype=np.float64 )
    W = np.asarray(W)
    W = W / W.sum(axis=0) # normalization step
    n = np.shape(W)[1]

    ## Create PPR matrix
    return beta*inv(np.eye(n)-(1.-beta)*W)

# Remove self-loops, multi-edges, and restrict to the largest component
def largest_component(G):
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from( selfLoops )
    return G.subgraph( sorted(nx.connected_components( G ), key=lambda cc: len(cc), reverse=True)[0] )
