#!/usr/bin/env python

# Load required modules
import networkx as nx, numpy as np, h5py
from scipy.linalg import inv, eigh
import hnio
from constants import *

################################################################################
# HOTNET2 DIFFUSION (PAGE RANK)
################################################################################

# Run the HotNet2 diffusion process on a given network
def hotnet2_diffusion(G, nodes, beta, verbose):
    if verbose > 1: print "* Creating HotNet2 diffusion matrix for beta=%s..." % beta
    W = nx.to_numpy_matrix( G , nodelist=nodes, dtype=np.float64 )
    W = np.asarray(W)
    W = W / W.sum(axis=0) # normalization step
    n = np.shape(W)[1]

    ## Create PPR matrix
    return beta*inv(np.eye(n)-(1.-beta)*W)

################################################################################
# HOTNET DIFFUSION (HEAT EQUATION)
################################################################################

# Run the HotNet diffusion process on a given network
def hotnet_diffusion(G, nodes, time, verbose):
    if verbose > 1: print "* Creating HotNet diffusion matrix for time t=%s..." % time
    L = nx.laplacian_matrix(G)
    Li = expm_eig( -time * L.todense() )
    return Li

def expm_eig(A):
    """
    Compute the matrix exponential for a square, symmetric matrix.
    """
    D, V = eigh(A)
    return np.dot(np.exp(D) * V, inv(V))

################################################################################
# HELPERS
################################################################################
        
# Run the entire HotNet2 diffusion process from start to finish
def save_diffusion_to_file( diffusion_type, diffusion_param, index_file, edge_file,
                            output_file, params=dict(), verbose=0):
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
    output = dict(edges=G.edges(), nodes=nodes)
    output.update(params.items())
    if diffusion_type == HOTNET2:
        output['beta'] = diffusion_param
        output['PPR']  = hotnet2_diffusion(G, nodes, output['beta'], verbose)
    elif diffusion_type == HOTNET:
        output['time'] = diffusion_param
        output['Li']   = hotnet_diffusion(G, nodes, output['time'], verbose)
    else:
        raise NotImplementedError('Diffusion of type "%s" not implemented' % diffusion_type)
        
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

# Remove self-loops, multi-edges, and restrict to the largest component
def largest_component(G):
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from( selfLoops )
    return G.subgraph( sorted(nx.connected_components( G ), key=lambda cc: len(cc), reverse=True)[0] )
