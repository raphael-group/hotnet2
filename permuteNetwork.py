#!/usr/bin/python

# Load required modules
import sys, networkx as nx
from hotnet2 import hnap

# Parse arguments
def parse_args(raw_args):
    description = 'Creates permuted versions of the given network, where each '\
                  'node retains its degree.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    # We set Q to be 115 by default (instead of the standard 100) because some
    # proportion of swaps fail as we are doing *connected* edge swaps
    parser.add_argument('-q', '--Q', default=115, type=float,
                        help='Edge swap constant.')
    parser.add_argument('-s', '--start_index', default=1, type=int,
                        help='Index to start output of permutations at.')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Edgelist file path.')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='Output prefix.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory.')
    parser.add_argument('-n', '--num_permutations', default=100, type=int,
                        help='Number of permuted networks to create.')

    return parser.parse_args(raw_args)

def permute_network( G, Q, num_edges, output_file ):
    # Permutes network by swapping edges Q * num_edges times
    H = G.copy()
    nswap = Q*num_edges
    swaps = nx.connected_double_edge_swap(H, nswap=nswap)
    nx.write_edgelist(H, output_file)
    return swaps

def run(args):
    # Load graph
    print "* Loading edge list.."
    G = nx.Graph()
    G.add_edges_from([ l.rstrip().split()[:2] for l in open(args.edgelist_file) ])
    num_edges, num_nodes = len(G.edges()), len(G.nodes())

    # Report info about the graph and number of swaps
    Q = args.Q
    print "\t- {} edges among {} nodes.".format(num_edges, num_nodes)
    print "\t- No. swaps to attempt = {}".format(Q*num_edges)

    # Permute graph the prescribed number of times
    print "* Creating permuted networks..."
    prefix, output_dir = args.output_prefix, args.output_dir
    swaps_performed = []
    increment = args.num_permutations/80. if args.num_permutations > 80 else 1
    for i in range(args.start_index, args.start_index + args.num_permutations):
        # Create a simple progress bar
        if int((i - args.start_index + 1) % increment) == 0:
            sys.stdout.write("+")
            sys.stdout.flush()

        # Create a new permuted network
        output_file = "{}/{}_edgelist_{}".format(output_dir, prefix, i)
        swaps = permute_network( G, Q, num_edges, output_file )
        swaps_performed.append( swaps * 1. )

    print

    # Report how many swaps were actually made
    avg_swaps = int(sum(swaps_performed) / float(len(swaps_performed)))
    print "* Avg. No. Swaps Made: {}".format(avg_swaps)
            
if __name__ == "__main__":
    run(parse_args(sys.argv[1:]))
