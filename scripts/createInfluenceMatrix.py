#!/usr/bin/python

# Load required modules
import sys, os, numpy as np, networkx as nx, scipy as sp, scipy.io
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
from hotnet2 import *

# Parse arguments
def get_parser():
    # Common arguments
    description = 'Create the personalized PageRank matrix for the given '\
                  'network and restart probability beta.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Location of edgelist file.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Location of gene-index file.')
    parser.add_argument('-n', '--network_name', required=True,
                        help='Name of network.')
    parser.add_argument('-o', '--output_file', required=True,
                    help="Output file.")
    parser.add_argument('-v', '--verbose', required=False, default=0, type=int, choices=range(5),
                    help="Control verbosity of output.")

    # Subparsers for the different diffusion types
    subparser = parser.add_subparsers(dest='diffusion_type', help='Diffusion method')

    hn2_parser = subparser.add_parser(HOTNET2)
    hn2_parser.add_argument('-b', '--beta', required=True, type=float, help="Restart probability beta.")

    hn_parser = subparser.add_parser(HOTNET)
    hn_parser.add_argument('-t', '--time', required=True, type=float, help='Diffusion time.')

    return parser

def run(args):
    params = dict(network_name=args.network_name)
    if args.diffusion_type == HOTNET2:
        diffusion_param = args.beta
    elif args.diffusion_type == HOTNET:
        diffusion_param = args.time
    else:
        raise NotImplementedError('Diffusion of type "%s" not implemented' % args.diffusion_type)
    save_diffusion_to_file(args.diffusion_type, diffusion_param, args.gene_index_file, args.edgelist_file,
                           args.output_file, params=params, verbose=args.verbose)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
