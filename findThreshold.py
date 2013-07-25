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
    
    #create parent parser for arguments common to both permutation types
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parent_parser.add_argument('-mn', '--infmat_name', default='Li',
                                help='Variable name of the influence matrices in the .mat files')
    parent_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                                help='Gene-index file for the influence matrices.')
    parent_parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    parent_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                                help='Number of permuted networks to use')
    parent_parser.add_argument('-s', '--test_statistic', choices=['max_cc_sizes', 'num_ccs'],
                               default='max_cc_sizes',
                               help='If max_cc_sizes, select smallest delta such that the size of\
                               the largest CC size is >= k. If num_ccs, select SMALLEST? delta that\
                               maximizes the number of CCs of size >= k.')
    parent_parser.add_argument('-l', '--sizes', nargs='+', type=int, help='See test_statistic')
    parent_parser.add_argument('-p', '--parallel', default=False, action='store_true',
                        help='Include flag to run permutation tests in parallel.')
    parent_parser.add_argument('-c', '--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    parent_parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    
    subparsers = parser.add_subparsers(title='Permutation techniques')

    #create subparser for options for permuting networks
    network_parser = subparsers.add_parser('network', help='Permute networks', parents=[parent_parser])
    network_parser.add_argument('-pnp', '--permuted_networks_path', required=True,
                                help='Path to influence matrices for permuted networks.\
                                      Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                                      path to be replaced with the iteration number')
    network_parser.set_defaults(delta_fn=get_deltas_for_network)
    
    #create subparser for options for permuting heat scores
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-mf', '--infmat_file', required=True,
                             help='Path to .mat file containing influence matrix')
    heat_parser.set_defaults(delta_fn=get_deltas_for_heat)
    
    #if l not specified, set default based on test statistic 
    args = parser.parse_args(raw_args)
    if not args.sizes:
        args.sizes = [5,10,15,20] if args.test_statistic == "max_cc_sizes" else [3,4,5]
                        
    return args

def run(args):
    heat, heat_params = hnio.load_heat_json(args.heat_file)
    deltas = args.delta_fn(args, heat)
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    args.delta_fn = args.delta_fn.__name__
    json.dump({"parameters": vars(args), "heat_parameters": heat_params,
               "deltas": deltas}, output_file, indent=4)
    if (args.output_file): output_file.close()


def get_deltas_for_network(args, heat):
    #construct list of paths to the first num_permutations     
    permuted_network_paths = [args.permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                              for i in range(1, args.num_permutations+1)]

    index2gene = hnio.load_index(args.infmat_index_file)

    deltas = delta.network_delta_selection(permuted_network_paths, args.infmat_name, index2gene,
                                           heat, args.sizes, not args.classic,
                                           args.parallel)
    
    return deltas


def get_deltas_for_heat(args, gene2heat):
    import scipy.io
    
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
    index2gene = hnio.load_index(args.infmat_index_file)
  
    M, gene_index = hn.induce_infmat(infmat, index2gene, sorted(gene2heat.keys()))

    heat_permutations = permutations.permute_heat(gene2heat, args.num_permutations,
                                                  parallel=args.parallel)
    delta_selection_fn = delta.find_best_delta_by_largest_cc if args.test_statistic == "max_cc_sizes" else delta.find_best_delta_by_num_ccs 
    
    deltas = delta.heat_delta_selection(M, gene_index, heat_permutations, args.sizes,
                                        not args.classic, args.parallel, delta_selection_fn)
    return deltas


if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
