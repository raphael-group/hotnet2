#!/usr/bin/env python

# Load required modules
import os, sys, argparse, json, networkx as nx
from itertools import combinations
from collections import defaultdict

# Argument parser
def get_parser():
    description = 'Constructs consensus subnetworks from HotNet(2) results.'

    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-d', '--directories', nargs='*', default=[],
        help='Paths to HotNet(2) results directories; the consensus procedure chooses from these directories.')
    parser.add_argument('-f', '--files', nargs='*', default=[],
        help='Paths to HotNet(2) results files.')
    parser.add_argument('-n', '--networks', nargs='*', required=True,
        help='Networks for HotNet(2) results directories or files.')
    parser.add_argument('-p', '--p_value_threshold', type=float, default=0.01,
        help='Threshold for p-values; default is 0.01.')
    parser.add_argument('-m', '--min_cc_size', type=int, default=2,
        help='Minimum connected component size; default is 2.')
    parser.add_argument('-o', '--output_file', required=True,
        help='Output file; provide .json extension for JSON output.')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Verbose')

    return parser

# Choose the results files when given directories of results files.
def choose_result_files(directories, min_cc_size, p_value_threshold, verbose=False):
    result_files = []
    # Consider each run.
    for directory in directories:
        result_statistics = []
        # Consider each \delta value for each run.
        for subdirectory in os.listdir(directory):
            if os.path.isdir(os.path.join(directory, subdirectory)) and subdirectory.startswith('delta'):
                results_file = os.path.join(directory, subdirectory, 'results.json')
                # For each \delt value, load the results for each \delta value and record the
                # number of component sizes k with p-values below a given threshold when k \geq
                # min_cc_size.
                with open(results_file, 'r') as f:
                    data = json.load(f)
                count = 0
                for k in data['statistics']:
                    if data['statistics'][k]['pval']>p_value_threshold and int(k)>=min_cc_size:
                        count += 1
                delta = data['parameters']['delta']
                result_statistics.append((count, delta, results_file))
        # Find the smallest \delta value with the largest number of component sizes k with p-values
        # below a threshold. For convenience, we find the smallest number of \delta values greater
        # greater than or equal to the threshold, and we sort from smallest to largest.
        selected_result_statistics = sorted(result_statistics)[0]
        if verbose:
            print '\t- Selected delta = {} for {}...'.format(selected_result_statistics[1], directory)
        result_files.append(selected_result_statistics[2])
    return result_files

# Load results.
def load_results(results_files, min_size):
    results = []
    for results_file in results_files:
        with open(results_file, 'r') as f:
            data = json.load(f)
        components = [ cc for cc in data['components'] if len(cc)>=min_size ]
        results.append(components)
    return results

# Construct consensus graph.
def consensus_edges(components, networks):
    edges_to_networks = defaultdict(set)
    for ccs, network in zip(components, networks):
        for cc in ccs:
            for u, v in combinations(cc, 2):
                edges_to_networks[(u, v)].add(network)
    return dict( (edge, len(edge_networks)) for edge, edge_networks in edges_to_networks.iteritems() )

def run(args):
    # Check arguments.
    if not args.directories and not args.files:
        raise ValueError('Neither the directories or files argument is specified; specify one argument.')
    if args.directories and args.files:
        raise ValueError('Both the directories and files arguments are specified; specify only one argument.')
    if args.directories and len(args.directories)!=len(args.networks):
        raise ValueError('The directories and networks arguments have different numbers of entries; provide equal numbers of entries.')
    if args.files and len(args.files)!=len(args.networks):
        raise ValueError('The files and networks arguments have different numbers of entries; provide equal numbers of entries.')

    # Count networks.
    num_networks = len(set(args.networks))
    if args.verbose:
        print '* Combining {} networks from {} HotNet(2) runs...'.format(num_networks, len(args.networks))

    # Choose results files if given directories of results files.
    if args.directories:
        if args.verbose:
            print '* Automatically choosing results from each run...'
        result_files = choose_result_files(args.directories, args.min_cc_size, args.p_value_threshold, args.verbose)
    else:
        if args.verbose:
            print '* Using provided results...'
        result_files = args.files

    # Load results.
    results = load_results(result_files, args.min_cc_size)

    # Create the full consensus graph.
    if args.verbose:
        print '* Constructing HotNet(2) consensus network...'
    edges = consensus_edges(results, args.networks)
    G = nx.Graph()
    G.add_weighted_edges_from( (u, v, w) for (u, v), w in edges.iteritems() )

    # Extract the connected components when restricted to edges in all networks.
    H = nx.Graph()
    H.add_edges_from( (u, v) for (u, v), w in edges.iteritems() if w >= num_networks )
    consensus = [ set(cc) for cc in nx.connected_components( H ) ]
    consensus_genes = set( g for cc in consensus for g in cc )

    # Expand each consensus by adding back any edges not in all networks.
    expanded_consensus = []
    linkers = set()
    for cc in consensus:
        other_consensus_genes = consensus_genes - cc
        neighbors = set( v for u in cc for v in G.neighbors(u) if v not in consensus_genes )
        expansion = set()
        for u in neighbors:
            cc_networks = max( G[u][v]['weight'] for v in set(G.neighbors(u)) & cc )
            consensus_neighbors = set( v for v in G.neighbors(u) if v in other_consensus_genes and G[u][v]['weight'] >= cc_networks )
            if len(consensus_neighbors) > 0:
                if any([ G[u][v]['weight'] > 1 for v in consensus_neighbors ]):
                    linkers.add( u )
            else:
                expansion.add( u )
        expanded_consensus.append( dict(consensus=list(cc), expansion=list(expansion)) )

    consensus_genes = set( g for cc in expanded_consensus for g in cc['consensus'] + cc['expansion'] )
    linkers -= consensus_genes

    # Summarize the results.
    if args.verbose:
        total_genes = set(v for ccs in results for cc in ccs for v in cc)
        print '* Returning {} consensus genes from {} original genes...'.format(len(consensus_genes), len(total_genes))

    # Output to file (either JSON or text depending on the file extension).
    with open(args.output_file, 'w') as f:
        if args.output_file.lower().endswith('.json'):
            # Convert the consensus to lists
            consensus = [ dict(core=list(c['consensus']), expansion=list(c['expansion'])) for c in expanded_consensus ]
            json.dump(dict(linkers=list(linkers), consensus=consensus), f, sort_keys=True, indent=4)
        else:
            output  = [ '# Linkers: {}'.format(', '.join(sorted(linkers))), '# Consensus:' ]
            output += [ '{}\t[{}] {}'.format(i, ', '.join(sorted(c['consensus'])), ', '.join(sorted(c['expansion']))) for i, c in enumerate(expanded_consensus) ]
            f.write( '\n'.join(output) )

if __name__ == '__main__':
    run(get_parser().parse_args( sys.argv[1:]) )
