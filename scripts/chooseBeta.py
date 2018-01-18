#!/usr/bin/python

# Load modules.
import numpy as np, scipy as sp, scipy.optimize
import sys, argparse

# Parse arguments.
def get_parser():
    description = 'Choose beta to balance the distribution of random walkers between neighbors and non-neighbors.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-d', '--digits', type=int, required=False, default=2, help='Number of digits in beta')
    parser.add_argument('-t', '--threshold', type=float, required=False, default=1.0, help='Threshold for edge weights')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Define functions; import these functions in package.
def load_edge_list(filename):
    '''
    Load edge list.
    '''
    edge_list = list()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()
                i = int(arrs[0])
                j = int(arrs[1])
                edge_list.append((i, j))

    return edge_list

def convert_edge_list_to_adjacency_matrix(edge_list):
    '''
    Convert an edge list to an adjacency matrix for an undirected graph.
    '''
    k = min(min(edge) for edge in edge_list)
    l = max(max(edge) for edge in edge_list)

    A = np.zeros((l-k+1, l-k+1), dtype=np.int)
    for i, j in edge_list:
        A[i-k, j-k] = A[j-k, i-k] = 1

    return A

def walk_matrix(A):
    '''
    Find the walk matrix for the random walk.
    '''
    d = np.sum(A, axis=0)
    d[np.where(d<=0)] = 1
    return np.asarray(A, dtype=np.float64)/np.asarray(d, dtype=np.float64)

def hotnet2_similarity_matrix(A, beta):
    '''
    Perform the random walk with restart process in HotNet2.
    '''
    from scipy.linalg import inv
    return beta*inv(np.eye(*np.shape(A))-(1-beta)*walk_matrix(A))

# Define more functions.
def difference(A, beta, threshold):
    '''
    Find difference between fraction of distribution on neighbors and non-neighbors.
    '''
    P = hotnet2_similarity_matrix(A, beta)
    np.fill_diagonal(P, 0)

    r = np.sum(P[np.where(A>=threshold)])
    s = np.sum(P[np.where(A<threshold)])
    return r-s

def balanced_beta(A, threshold, digits):
    '''
    Find value of beta that sets difference to zero between fraction of distribution on neighbors and non-neighbors to zero.
    '''
    return sp.optimize.ridder(lambda beta: difference(A, beta, threshold), a=0.1**digits, b=1.0-0.1**digits, xtol=0.1**(digits+1))

# Run script.
def run(args):
    # Load graph.
    edge_list = load_edge_list(args.edge_list_file)
    A = convert_edge_list_to_adjacency_matrix(edge_list)

    # Compute beta.
    beta = balanced_beta(A, args.threshold, args.digits)

    # Save beta.
    fmt = '{:.' + str(args.digits) + 'f}'
    with open(args.output_file, 'w') as f:
        f.write(fmt.format(beta))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
