# -*- coding: iso-8859-1 -*-
import hnio
import hotnet2 as hn
import delta
import permutations
import networkx as nx
import sys
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
    #TODO fix this so that multithreading default is true
    parser.add_argument('-m', '--multithreaded', default=False, action='store_true',
                        help='Set to 0 to disable running permutation tests in parallel')
    parser.add_argument('--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')

    subparsers = parser.add_subparsers(title='Permutation techinques')

    #create subparser for options for permuting networks
    network_parser = subparsers.add_parser('network', description='Permute networks')
    network_parser.add_argument('-p', '--permuted_networks_path', required=True,
                                help='Path to influence matrices for permuted networks.\
                                      Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                                      path to be replaced with the iteration number')
    network_parser.add_argument('-mn', '--infmat_name', default='Li',
                                help='Variable name of the influence matri[x/ces] in the .mat file[s]')
    network_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                                help='Gene-index file for the influence matrix.')
    network_parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    network_parser.add_argument('-k', '--max_cc_sizes', nargs='+', type=int, required=True, 
                                help='Max CC sizes for delta selection')
    network_parser.add_argument('-n', '--num_permutations', type=int,
                                help='Number of permuted networks to use')
    network_parser.set_defaults(func=run_for_network)
    
    #create subparser for options for permuting heat scores
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores')
    heat_parser.add_argument('-mf', '--infmat_file', required=True,
                             help='Path to .mat file containing influence matrix')
    heat_parser.add_argument('-mn', '--infmat_name', required=True, default='Li',
                             help='Variable name of the influence matrix in the .mat file')
    heat_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                             help='Gene-index file for the influence matrix.')
    heat_parser.add_argument('-hf', '--heat_file', required=True,
                             help='Heat score file')
    #TODO: make k and l mutually exclusive
    heat_parser.add_argument('-l', '--max_cc_sizes', nargs='+', type=int, required=True, 
                             help='Max CC sizes for delta selection')
    # heat_parser.add_argument('-k', '--test_cc_size', nargs='+', type=int, required=True, 
    #                          help='Value for choosing delta to maximize # CCs of size >= k')
    heat_parser.add_argument('-n', '--num_permutations', type=int,
                             help='Number of heat score permutations to test')
    heat_parser.set_defaults(func=run_for_heat)
    #TODO add ability to specify pre-permuted heat files <- is this actually useful?
                        
    return parser.parse_args(raw_args)


def run_for_network(args):
    #construct list of paths to the first num_permutations     
    permuted_network_paths = [args.permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i)) for i in range(1, args.num_permutations+1)]

    index2gene = hnio.load_index(args.infmat_index_file)
    heat = hnio.load_heat(args.heat_file)

    h_vec = [heat[gene] for index, gene in sorted(index2gene.items()) if gene in heat]
    component_fn = strong_ccs if not args.classic else nx.connected_components
        
    #TODO: at some point, pass around immutable views -- deferring for now since there's no built-in immutable dict type
    deltas = delta.network_delta_selection(permuted_network_paths, index2gene, args.infmat_name, sorted(heat.keys()),
                                           h_vec, args.max_cc_sizes, component_fn, args.multithreaded)
    
    print "Deltas is: ", deltas


def run_for_heat(args):
    import scipy.io
    
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
    index2gene = hnio.load_index(args.infmat_index_file)
    gene2heat = hnio.load_heat(args.heat_file)
  
    M, gene_index = hn.induce_infmat(infmat, index2gene, sorted(gene2heat.keys()))

    heat_permutations = permutations.permute_heat(gene2heat, args.num_permutations, parallel=args.multithreaded)
    deltas = delta.heat_delta_selection(M, gene_index, heat_permutations, args.max_cc_sizes, not args.classic, args.multithreaded)
    print deltas


if __name__ == "__main__": 
    args = parse_args(sys.argv[1:])
    args.func(args)
