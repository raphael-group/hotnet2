#!/usr/bin/python

# Load required modules
import sys, os, numpy as np, networkx as nx, scipy as sp, scipy.io
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
from hotnet2 import * 

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
    parser.add_argument('-b', '--beta', required=True, type=float,
                    help="Restart probability beta.")
    parser.add_argument('--exclude_network', required=False, action='store_true',
                    default=False,
                    help="Exclude the gene index and edges in the output file.")
    parser.add_argument('--verbose', required=False, default=0, type=int, choices=range(5),
                    help="Control verbosity of output.")
    return parser

def run(args):
    params = dict(network_name=args.network_name, beta=beta, permuted_networks_path=args.permuted_networks_path)
    save_hotnet2_diffusion_to_file(args.gene_index_file, args.edge_list_file, args.beta, args.output_file, args.verbose)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
