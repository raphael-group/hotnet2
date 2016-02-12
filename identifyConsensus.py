#!/usr/bin/env python

# Load required modules
import sys, json, argparse, networkx as nx, os
from collections import defaultdict

# Argument parser
def get_parser():
    description = 'Constructs consensus subnetworks from HotNet(2) results.'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-r', '--results_files',
                        nargs="*", help='Paths to HotNet(2) results files (/directories if -d flag is set).', required=True)
    parser.add_argument('-n', '--networks', nargs="*", required=True,
                        help='List the network that corresponds to each results file (/directory if -d flag is set).')
    parser.add_argument('-d', '--delta_choice', action='store_true',
                        help='Automated delta selection. If flag is set, -r will expect a single directory per HotNet2 run.')
    parser.add_argument('-o', '--output_file',required=True,
                        help='Output file. Text output by default. Use a .json extension to get JSON output.')
    parser.add_argument('-ms', '--min_cc_size', help='Min CC size.', type=int, default=2)

    return parser

# Choose results files to use for consensus (1 per network)
def choose_results( results_files, networks, min_cc_size, delta_choice ):
    if delta_choice:
        networks_results = defaultdict( lambda: set(), dict() )
        network_dirs = results_files
        results_files = []
        for network_dir, network in zip( network_dirs, networks ):
            results_dirs = [ os.path.join(network_dir, dir)
                             for dir in os.listdir(network_dir)
                             if os.path.isdir(os.path.join(network_dir, dir)) & ( dir.split('_')[0] == 'delta' ) ]
            for results_dir in results_dirs:
                networks_results[network].add( results_dir + "/results.json" )
        for network in networks:
            delta_list = []
            # Loading results json files and adding entries if k (size) > min_cc_size into a list of
            # tuples: ( delta, k, p-value, results json file )
            for results_file in networks_results[network]:
                with open( results_file ) as f: obj = json.load( f )
                result_delta_list = [ ( float( obj['parameters']['delta'] ), int( k ), stats['pval'], results_file )
                                      for k, stats in obj['statistics'].iteritems()
                                      if int( k ) >= min_cc_size ]
                delta_list.extend( result_delta_list )
            # Choosing results json file with the lowest p-value, then the lowest k, then the lowest delta
            results_files.append( sorted( delta_list )[0][3] )

    return results_files, networks

def load_results(results_files, min_cc_size):
    results = []
    for filename in results_files:
        # Load the components from each results file
        print '* Loading {}...'.format(filename)
        with open(filename) as f: obj = json.load(f)
        components = [cc for cc in obj['components'] if len(cc) >= min_cc_size]
        results.append( components )

    return results

# Construct consensus graph
def consensus_edges(results, networks, all_genes):
    # Create membership dictionary
    gene2neighbors = dict((g, defaultdict(set)) for g in all_genes)
    for ccs, network in zip(results, networks):
        for cc in ccs:
            for i, u in enumerate(cc):
                for j, v in enumerate(cc):
                    if i != j:
                        gene2neighbors[u][v].add( network )

    # Create an edge list from the membership dictionary
    from itertools import combinations
    edges = [(g1, g2, dict(networks=len(gene2neighbors[g1][g2])))
             for g1, g2 in combinations(all_genes, 2)
             if len(gene2neighbors[g1][g2]) > 0 ]

    return edges

def run(args):
    # Assert that a network was given for each results file
    if len(args.results_files) != len(args.networks):
        raise ValueError("You must pass in one network name for each result file.")

    num_networks = len(set(args.networks))

    # Choose results to use
    results_files, networks = choose_results( args.results_files, args.networks, args.min_cc_size, args.delta_choice )

    # Load results
    print "* Loading results..."
    results = load_results(results_files, args.min_cc_size)

    # Create a set that lists all genes found in the HotNet2 results
    all_genes = set( g for ccs in results for cc in ccs for g in cc )
    print "\t- %s genes included across all results" % len(all_genes)

    # Create the full consensus graph
    edges = consensus_edges( results, networks, all_genes )
    G = nx.Graph()
    G.add_edges_from( edges )

    # Extract the connected components when restricted to edges in all networks
    H = nx.Graph()
    H.add_edges_from( (u, v, d) for u, v, d in edges if d['networks'] >= num_networks )
    consensus = [ set(cc) for cc in nx.connected_components( H ) ]
    consensus_genes = set( g for cc in consensus for g in cc )

    # Expand each consensus by adding back any edges not in all networks
    expanded_consensus = []
    linkers = set()
    for cc in consensus:
        other_consensus_genes = consensus_genes - cc
        neighbors = set( v for u in cc for v in G.neighbors(u) if v not in consensus_genes )
        expansion = set()
        for u in neighbors:
            cc_networks = max( G[u][v]['networks'] for v in set(G.neighbors(u)) & cc )
            consensus_neighbors = set( v for v in G.neighbors(u) if v in other_consensus_genes and G[u][v]['networks'] >= cc_networks )
            if len(consensus_neighbors) > 0:
                if any([ G[u][v]['networks'] > 1 for v in consensus_neighbors ]):
                    linkers.add( u )
            else:
                expansion.add( u )
        expanded_consensus.append( dict(consensus=list(cc), expansion=list(expansion)) )

    consensus_genes = set( g for cc in expanded_consensus for g in cc['consensus'] + cc['expansion'] )
    linkers -= consensus_genes

    print "* No. consensus genes:", len(consensus_genes)

    # Output to file (JSON or text depending on the file extension)
    with open(args.output_file, "w") as out:
        if args.output_file.lower().endswith(".json"):
            # Convert the consensus to lists
            consensus = [ dict(core=list(c['consensus']), expansion=list(c['expansion'])) for c in expanded_consensus ]
            json.dump(dict(linkers=list(linkers), consensus=consensus), out, sort_keys=True, indent=4)
        else:
            output  = [ "# Linkers: {}".format(", ".join(sorted(linkers))), "#Consensus" ]
            output += ["{}\t[{}] {}".format(i, ", ".join(sorted(c['consensus'])), ", ".join(sorted(c['expansion']))) for i, c in enumerate(expanded_consensus) ]
            out.write( "\n".join(output) )

if __name__ == "__main__": run(get_parser().parse_args(sys.argv[1:]))
