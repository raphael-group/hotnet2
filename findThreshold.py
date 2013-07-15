# -*- coding: iso-8859-1 -*-
import hnap
import hnio
import hotnet2 as hn
import delta
import permutations
import sys
import json

ITERATION_REPLACEMENT_TOKEN = '##NUM##'

def parse_args(raw_args):
    description = "Runs HotNet threshold-finding procedure.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python findThreshold.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-p', '--parallel', default=False, action='store_true',
                        help='Include flag to run permutation tests in parallel.')
    parser.add_argument('-c', '--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')

    subparsers = parser.add_subparsers(title='Permutation techinques')

    #create subparser for options for permuting networks
    network_parser = subparsers.add_parser('network', help='Permute networks')
    network_parser.add_argument('-pnp', '--permuted_networks_path', required=True,
                                help='Path to influence matrices for permuted networks.\
                                      Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                                      path to be replaced with the iteration number')
    network_parser.add_argument('-mn', '--infmat_name', default='Li',
                                help='Variable name of the influence matrices in the .mat files')
    network_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                                help='Gene-index file for the influence matrices.')
    network_parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    network_parser.add_argument('-k', '--max_cc_sizes', nargs='+', type=int, default=[5,10,15,20], 
                                help='Max CC sizes for delta selection')
    network_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                                help='Number of permuted networks to use')
    network_parser.set_defaults(delta_fn=get_deltas_for_network)
    
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
    heat_parser.add_argument('-l', '--max_cc_sizes', nargs='+', type=int, default=[5,10,15,20],
                             help='Max CC sizes for delta selection')
    #TODO: make k and l mutually exclusive
    # heat_parser.add_argument('-k', '--test_cc_size', nargs='+', type=int, required=True, 
    #                          help='Value for choosing delta to maximize # CCs of size >= k')
    heat_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                             help='Number of heat score permutations to test')
    heat_parser.set_defaults(delta_fn=get_deltas_for_heat)
                        
    return parser.parse_args(raw_args)

def run(args):
    deltas = args.delta_fn(args)
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    args.delta_fn = args.delta_fn.__name__
    json.dump({"parameters": vars(args), "deltas": deltas}, output_file, indent=4)
    if (args.output_file): output_file.close()


def get_deltas_for_network(args):
    #construct list of paths to the first num_permutations     
    permuted_network_paths = [args.permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                              for i in range(1, args.num_permutations+1)]

    index2gene = hnio.load_index(args.infmat_index_file)
    heat = hnio.load_heat(args.heat_file)

    deltas = delta.network_delta_selection(permuted_network_paths, args.infmat_name, index2gene,
                                           heat, args.max_cc_sizes, not args.classic,
                                           args.parallel)
    
    return deltas


def get_deltas_for_heat(args):
    import scipy.io
    
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
    index2gene = hnio.load_index(args.infmat_index_file)
    gene2heat = hnio.load_heat(args.heat_file)
  
    M, gene_index = hn.induce_infmat(infmat, index2gene, sorted(gene2heat.keys()))

    heat_permutations = permutations.permute_heat(gene2heat, args.num_permutations,
                                                  parallel=args.parallel)
    deltas = delta.heat_delta_selection(M, gene_index, heat_permutations, args.max_cc_sizes,
                                        not args.classic, args.parallel)
    return deltas


if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
