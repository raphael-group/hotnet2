# -*- coding: iso-8859-1 -*-
import os.path
from sys import argv
from delta import *
import networkx as nx
strong_ccs = nx.strongly_connected_components

ITERATION_REPLACEMENT_TOKEN = '##NUM##'

def parse_args(raw_args):  
    import argparse
    
    class BetterFileArgParser(argparse.ArgumentParser):
        def convert_arg_line_to_args(self, arg_line):
            if not arg_line.startswith('#'):
                for arg in arg_line.split():
                    yield arg
            
    description = "" #TODO
    parser = BetterFileArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-p', '--permuted_networks_path', required=True,
                        help='Path to influence matrices for permuted networks.\
                              Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                              path to be replaced with the iteration number')
    parser.add_argument('-mn', '--infmat_name', default='Li',
                        help='Variable name of the influence matrices in the .mat files')
    parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                        help='Gene-index file for the influence matrix.')
    parser.add_argument('-n', '--num_permutations', type=int,
                        help='Number of permuted networks to use')
    parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    parser.add_argument('-k', '--max_cc_sizes', nargs='+', type=int, required=True, 
                        help='Max CC sizes for delta selection')
    #TODO fix this so that multithreading default is true
    parser.add_argument('-m', '--multithreaded', default=False, action='store_true',
                        help='Set to 0 to disable running permutation tests in parallel')
    parser.add_argument('--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
                        
    return parser.parse_args(raw_args)

def run(args):
    #construct list of paths to the first num_permutations     
    permuted_network_paths = [args.permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i)) for i in range(1, args.num_permutations+1)]

    index2gene = load_index(args.infmat_index_file)
    heat = load_heat(args.heat_file)

    h_vec = [heat[gene] for index, gene in sorted(index2gene.items()) if gene in heat]
    component_fn = strong_ccs if not args.classic else nx.connected_components
        
    #TODO: at some point, pass around immutable views -- deferring for now since there's no built-in immutable dict type
    deltas = network_delta_selection(permuted_network_paths, index2gene, args.infmat_name, sorted(heat.keys()),
                                     h_vec, args.max_cc_sizes, component_fn, args.multithreaded)
    #def network_delta_selection(permuted_network_paths, index2gene, infmat_name, tested_genes, h, sizes, component_fn=strong_ccs, parallel=True):
    
    print "Deltas is: ", deltas
    #max_sizes = []
    #for delta in deltas[15]:
        #max_sizes.append( max( component_sizes( weighted_graph( sim, gene_index, delta
            #) ) ))
    #print "Max sizes is: ", max_sizes

    #open("delta_selection.tsv", 'w').write(
        #"\n".join(["%s\t%s" % (delta, size) for delta, size in zip(deltas,
            #max_sizes)])
            #)

if __name__ == "__main__": 
    run( parse_args(argv[1:]) )