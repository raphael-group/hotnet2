#!/usr/bin/python

# Import required modules
import sys, os, json, networkx as nx
from collections import defaultdict

def parse_args(input_list=None):
    # Parse args
    import argparse
    class Args: pass
    args = Args()
    description = 'Constructs consensus subnetworks from HotNet(2) results.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--result_directories', '-n',
                        nargs="*", help='Networks.', required=True)
    parser.add_argument('-o', '--output_file', help='Output file.', required=True)
    parser.add_argument('-ms', '--min_cc_size', help='Min CC size.', type=int, default=2)

    if not input_list: parser.parse_args(namespace=args)
    else: parser.parse_args(input_list, namespace=args)

    return args
    
def load_results(networks, min_cc_size):
    results = dict()
    for n in networks:
        results[n] = dict()
        print '* Loading {}...'.format(n)
        # Load the results
        results_file = "{}/results.json".format(n)
        obj = json.load(open(results_file))

        # Load the components
        components = [cc for cc in obj['components'] if len(cc) >= min_cc_size]

        # Store the results
        results[n] = components

    return results

# Construct consensus graph
def consensus_edges( results, all_genes ):
    # Create membership dictionary
    def new_neighbors(): return dict(networks=0)
    gene2neighbors = dict([(g, defaultdict(new_neighbors)) for g in all_genes])    
    for n in results:
        for cc in results[n]:
            for i in range(len(cc)):
                for j in range(len(cc)):
                    if i != j:
                        gene2neighbors[cc[i]][cc[j]]['networks'] += 1
                            
    # Create an edge list from the membership dictionary
    from itertools import combinations
    edges = [(g1, g2, dict(networks=gene2neighbors[g1][g2]['networks']))
             for g1, g2 in combinations(all_genes, 2)
             if gene2neighbors[g1][g2]['networks'] > 0 ]

    return edges

def run( args ):
    # Load the deltas
    networks = args.result_directories
    # Load results
    print "* Loading results..."
    results = load_results(networks, args.min_cc_size)

    # Create a set that lists all genes found in the HotNet2 results
    all_genes = set([g for n in results
                     for cc in results[n] for g in cc])
    print "\t- %s genes included across all results" % len(all_genes)
            
    # Create the full consensus graph
    edges = consensus_edges( results, all_genes )
    G = nx.Graph()
    G.add_edges_from( edges )

    # Extract the connected components when restricted to edges with 3 networks
    H = nx.Graph()
    H.add_edges_from( (u, v, d) for u, v, d in edges if d['networks'] == 3 )
    consensus = [ set(cc) for cc in nx.connected_components( H ) ]
    consensus_genes = set( g for cc in consensus for g in cc )
    
    # Expand each consensus by adding back any edges with < 3 networks
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

    # Output to file
    output  = [ "# Linkers: {}".format(", ".join(sorted(linkers))), "#Consensus" ]
    output += ["{}\t[{}] {}".format(i, ", ".join(sorted(c['consensus'])), ", ".join(sorted(c['expansion']))) for i, c in enumerate(expanded_consensus) ]
    open(args.output_file, "w").write( "\n".join(output) )

    # Annotate the subnetworks each linker came from
    for g in linkers:
        neighbors = set( u for u in G.neighbors(g) if G[g][u]['networks'] > 1 )
        linked = [ "[{}]".format(" ".join(cc['consensus'])) for i, cc in enumerate(expanded_consensus) if len(set(cc['consensus'] + cc['expansion']) & neighbors) != 0 ]
        print g, ", ".join(map(str, linked))
            
            

if __name__ == "__main__": run( parse_args() )
