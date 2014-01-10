# -*- coding: iso-8859-1 -*-
import sys
import json
import scipy.io
from hotnet import hnap, hnio, delta, permutations
from hotnet.constants import *

def parse_args(raw_args):
    description = "Runs HotNet threshold-finding procedure.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python findThreshold.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mn', '--infmat_name', default='PPR',
                        help='Variable name of the influence matrices in the .mat files')
    parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                        help='Path to tab-separated file containing an index in the first\
                              column and the name of the gene represented at that index in\
                              the second column of each line.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='JSON heat score file generated via generateHeat.py')
    parser.add_argument('-pnp', '--permuted_networks_path', required=True,
                        help='Path to influence matrices for permuted networks.\
                              Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                              path to be replaced with the iteration number')
    parser.add_argument('-n', '--num_permutations', type=int, required=True,
                        help='Number of permuted data sets to generate')
    parser.add_argument('-l', '--sizes', nargs='+', type=int, default=[5,10,15,20],
                        help='Largest CC size')
    parser.add_argument('--parallel', dest='parallel', action='store_true',
                        help='Run permutation tests in parallel. Only recommended for machines\
                              with at least 8 cores.')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                        help='Run permutation tests sequentially. Recommended for machines\
                              with fewer than 8 cores.')
    parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    parser.set_defaults(parallel=False)
    
    return parser.parse_args(raw_args)

def run(args):  
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)
    deltas = get_deltas_for_network(args.permuted_networks_path, heat, args.infmat_name,
                                    infmat_index, args.sizes, args.num_permutations, args.parallel)
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump({"parameters": vars(args), "heat_parameters": heat_params,
               "deltas": deltas}, output_file, indent=4)
    if (args.output_file): output_file.close()

def get_deltas_for_network(permuted_networks_path, heat, infmat_name, index2gene, sizes,
                           num_permutations, parallel):
    print "* Performing permuted network delta selection..."
    
    #construct list of paths to the first num_permutations permutations
    permuted_network_paths = [permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                              for i in range(1, num_permutations+1)]

    return delta.network_delta_selection(permuted_network_paths, infmat_name, index2gene, heat,
                                         sizes, parallel, delta.find_best_delta_by_largest_cc)

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
