#!/usr/bin/env python

# Load required modules
import os, sys, argparse, json, networkx as nx
from itertools import combinations
from collections import defaultdict

# Load HotNet2
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
from hotnet2.consensus import identify_consensus

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

# Find the HotNet2 results file in each directory, grouping each file
# by directory
def find_results_files(directories, verbose=0):
    # Identify the results files
    results_files = []
    for directory in directories:
        # Consider each \delta value for each run.
        result_group = []
        for subdirectory in os.listdir(directory):
            if os.path.isdir(os.path.join(directory, subdirectory)) and subdirectory.startswith('delta'):
                result_group.append(os.path.join(directory, subdirectory, 'results.json'))
        results_files.append(result_group)

    return results_files

# Choose the results files when given directories of results files.
def load_single_runs(results_groups, verbose=0):
    # Load each run
    single_runs = []
    for result_group in results_groups:
        # We are assuming all the results in the same directory will be
        # for the same heat score and network
        run_group = []
        for results_file in result_group:
            with open(results_file, 'r') as f:
                data = json.load(f)
                heat_name = data['parameters']['heat_name']
                network_name = data['parameters']['network_name']
                run_group.append( (data['components'], data['statistics'], data['parameters']['delta']) )
        single_runs.append((network_name, heat_name, run_group))

    return single_runs

def run(args):
    # Check arguments.
    if not args.directories and not args.files:
        raise ValueError('Neither the directories or files argument is specified; specify one argument.')
    if args.directories and args.files:
        raise ValueError('Both the directories and files arguments are specified; specify only one argument.')

    # Choose results files if given directories of results files.
    if args.directories:
        result_files = find_results_files(args.directories, args.verbose)
    else:
        result_files = [ [f] for f in args.files ]

    # Load results.
    single_runs = load_single_runs(result_files, args.verbose)

    # Create the full consensus graph.
    consensus, linkers, auto_deltas = identify_consensus(single_runs)
    consensus_genes = set( g for cc in consensus for g in cc['core'] + cc['expansion'] )

    # Summarize the results.
    if args.verbose:
        total_genes = set(v for ccs in results for cc in ccs for v in cc)
        print '* Returning {} consensus genes from {} original genes...'.format(len(consensus_genes), len(total_genes))

    # Output to file (either JSON or text depending on the file extension).
    with open(args.output_file, 'w') as f:
        if args.output_file.lower().endswith('.json'):
            # Convert the consensus to lists
            consensus = [ dict(core=list(c['core']), expansion=list(c['expansion'])) for c in consensus ]
            output    = dict(linkers=list(linkers), consensus=consensus, deltas=auto_deltas)
            json.dump(output, f, sort_keys=True, indent=4)
        else:
            output  = [ '# Linkers: {}'.format(', '.join(sorted(linkers))), '# Core\tExpansion' ]
            output += [ '{}\t{}'.format(', '.join(sorted(c['core'])), ', '.join(sorted(c['expansion']))) for c in consensus ]
            f.write( '\n'.join(output) )

if __name__ == '__main__':
    run(get_parser().parse_args( sys.argv[1:]) )
