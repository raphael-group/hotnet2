# -*- coding: iso-8859-1 -*-
from sys import argv
from stats import *
import json

def parse_args(raw_args):  
    import argparse
    
    class BetterFileArgParser(argparse.ArgumentParser):
        def convert_arg_line_to_args(self, arg_line):
            if not arg_line.startswith('#'):
                for arg in arg_line.split():
                    yield arg
            
    description = "" #TODO
    parser = BetterFileArgParser(description=description, fromfile_prefix_chars='@')

    #TODO: add option to give permutation files as input
    parser.add_argument('-hof', '--hotnet_output_file', required=True,
                        help='Path to output file produced by runHotnet2.py')
    parser.add_argument('-pgf', '--permutation_genes_file',
                        help='Path to file containing a list of additional genes that can have\
                              permuted heat values assigned to them in permutation tests')
    parser.add_argument('-n', '--num_permutations', type=int, required=True,
                        help='Number of permutation tests to run')
    parser.add_argument('-s', '--cc_start_size', type=int, default=2,
                        help='Smallest connected component size to count')
    parser.add_argument('-l', '--cc_stop_size', type=int, default=10,
                        help='Largest connected component size to count')
    #TODO: fix this so that multithreading default is true
    parser.add_argument('-m', '--multithreaded', default=False, action='store_true',
                        help='Set to 0 to disable running permutation tests in parallel')
                        
    return parser.parse_args(raw_args)

PARAMETERS = "parameters"

def run(args):
    import scipy.io

    hotnet_output = json.load(open(args.hotnet_output_file))

    infmat = scipy.io.loadmat(hotnet_output[PARAMETERS]["infmat_file"])[hotnet_output[PARAMETERS]["infmat_name"]]  
    infmat_index = load_index(hotnet_output[PARAMETERS]["infmat_index_file"])
  
    #heat is a dict from gene names to heat scores
    heat = load_heat(hotnet_output[PARAMETERS]["heat_file"])

    delta = hotnet_output[PARAMETERS]["delta"]
    sizes = range(args.cc_start_size, args.cc_stop_size+1)
  
    M, gene_index, inf_score = induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = heat_vec(heat, gene_index)
    sim, sim_score = similarity_matrix(M, h, gene_index, not hotnet_output[PARAMETERS]["classic"])
    G = weighted_graph(sim, gene_index, delta)

    print "* Sizes:"
    print component_sizes_str( G, 4 )
    print "* Components:"
    print write_components( G, 4 )

    extra_genes = load_permutation_genes(args.permutation_genes_file)

    genes_eligible_for_heat = sorted([g for g in (set(gene_index.values()) | extra_genes) if g in infmat_index.values()])
    print len(genes_eligible_for_heat)

     # size2stats = compute_significance( DiG, infmat, infmat_index, genes_in_network,
     #                                    h, DELTA, filtered_genes_in_network, 1, parallel=False )

#     #size2counts is dict(size -> list of counts, 1 per permutation)
    
    sizes2counts = calculate_permuted_cc_counts(infmat, infmat_index, genes_eligible_for_heat, h, delta,
                                                sorted(set(gene_index.values())), args.num_permutations,
                                                sizes, not hotnet_output[PARAMETERS]["classic"],
                                                args.multithreaded)


    #compute_significance( DiG, infmat, infmat_index, genes_in_network,
#                                    h, DELTA, filtered_genes_in_network, 1, parallel=False )
    real_counts = num_components_min_size(G, sizes)
    size2real_counts = dict(zip(sizes, real_counts))
    sizes2stats = compute_statistics(size2real_counts, sizes2counts, args.num_permutations)
    print json.dumps(sizes2stats, indent=4)

if __name__ == "__main__": 
    run( parse_args(argv[1:]) )