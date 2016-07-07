#!/usr/bin/python

# Load required modules
import sys, os, networkx as nx, multiprocessing as mp
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))

# Parse arguments
def get_parser():
    from hotnet2 import hnap
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
    parser.add_argument('-c', '--cores', default=1, type=int,
                        help='Use given number of cores. Pass -1 to use all available.')

    return parser

def permute_network( G, Q, numEdges, outputFile ):
    # Permutes network by swapping edges Q * numEdges times
    H = G.copy()
    nswap = Q*numEdges
    swaps = nx.connected_double_edge_swap(H, nswap=nswap)
    nx.write_edgelist(H, outputFile)
    return swaps

store = dict(maxSeen=-1)
def permute_network_wrapper((G, Q, numEdges, outputFile, i, n)):
    swaps = permute_network( G, Q, numEdges, outputFile )
    store['maxSeen'] = max(store['maxSeen'], i)
    sys.stdout.write("\r{}/{}".format(store['maxSeen'], n))
    sys.stdout.flush()
    return swaps

def run(args):
    # Load graph
    from hotnet2 import largest_component
    print "* Loading edge list.."
    G = nx.Graph()
    with open(args.edgelist_file) as infile:
        G.add_edges_from([ l.rstrip().split()[:2] for l in infile ])
        G = largest_component(G)
    numEdges, numNodes = len(G.edges()), len(G.nodes())

    # Report info about the graph and number of swaps
    Q = args.Q
    print "\t- {} edges among {} nodes.".format(numEdges, numNodes)
    print "\t- No. swaps to attempt = {}".format(Q*numEdges)

    # Permute graph the prescribed number of times
    print "* Creating permuted networks..."
    def outputFileName(i):
        return "{}/{}_edgelist_{}".format(args.output_dir, args.output_prefix, i+args.start_index)

    n = args.num_permutations
    if args.cores != 1:
        cores = mp.cpu_count() if args.cores == -1 else min(args.cores, mp.cpu_count)
        jobArgs = [ (G, Q, numEdges, outputFileName(i), i+args.start_index, n) for i in range(n) ]
        pool    = mp.Pool(cores)
        swaps = pool.map(permute_network_wrapper, jobArgs)
        pool.close()
        pool.join()
    else:
        swaps = [ permute_network_wrapper((G, Q, numEdges, outputFileName(i), i+args.start_index, n))
                  for i in range(n) ]

    print

    # Report how many swaps were actually made
    avgSwaps = int(sum(swaps) / float(len(swaps)))
    print "* Avg. No. Swaps Made: {}".format(avgSwaps)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
