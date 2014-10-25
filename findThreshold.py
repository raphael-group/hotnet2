# -*- coding: iso-8859-1 -*-
import json
import scipy.io
import sys
from hotnet2 import delta, hnap, hnio, permutations
from hotnet2.constants import MAX_CC_SIZE, NUM_CCS, ITERATION_REPLACEMENT_TOKEN

def get_parser():
    description = "Runs HotNet threshold-finding procedure.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python findThreshold.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    #create parent parser for arguments common to both permutation types
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parent_parser.add_argument('-mn', '--infmat_name', default='PPR',
                               help='Variable name of the influence matrices in the .mat files')
    parent_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                               help='Path to tab-separated file containing an index in the first\
                                     column and the name of the gene represented at that index in\
                                     the second column of each line.')
    parent_parser.add_argument('-hf', '--heat_file', required=True,
                               help='JSON heat score file generated via generateHeat.py')
    parent_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                                help='Number of permuted data sets to generate')
    parent_parser.add_argument('-s', '--test_statistic', choices=[MAX_CC_SIZE, NUM_CCS],
                               default=MAX_CC_SIZE,
                               help='If ' + MAX_CC_SIZE +', select smallest delta for each permuted\
                                     dataset such that the size of the largest CC is <= l. If ' +
                                     NUM_CCS + 'select for each permuted dataset the delta that \
                                     maximizes the number of CCs of size >= l.')
    parent_parser.add_argument('-l', '--sizes', nargs='+', type=int, help='See test_statistic')
    parent_parser.add_argument('--parallel', dest='parallel', action='store_true',
                               help='Run permutation tests in parallel.')
    parent_parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                               help='Run permutation tests sequentially.')
    parent_parser.add_argument('-c', '--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    parent_parser.add_argument('-o', '--output_file',
                               help='Output file.  If none given, output will be written to stdout.')
    parent_parser.set_defaults(parallel=False)
    
    subparsers = parser.add_subparsers(title='Permutation techniques', dest="perm_type")

    #create subparser for options for permuting networks
    network_parser = subparsers.add_parser('network', help='Permute networks', parents=[parent_parser])
    network_parser.add_argument('-pnp', '--permuted_networks_path', required=True,
                                help='Path to influence matrices for permuted networks.\
                                      Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                                      path to be replaced with the iteration number')
    
    #create subparser for options for permuting heat scores
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-mf', '--infmat_file', required=True,
                             help='Path to .mat file containing influence matrix')
    heat_parser.add_argument('-pgf', '--permutation_genes_file', default=None,
                             help='Path to file containing a list of additional genes that can have\
                                   permuted heat values assigned to them in permutation tests')
    
    #create subparser for options for permuting mutation data
    mutation_parser = subparsers.add_parser('mutations', help='Permute mutation data',
                                            parents=[parent_parser])
    mutation_parser.add_argument('-mf', '--infmat_file', required=True,
                                 help='Path to .mat file containing influence matrix')
    mutation_parser.add_argument('-glf', '--gene_length_file', required=True,
                                 help='Path to tab-separated file containing gene names in the\
                                       first column and the length of the gene in base pairs in\
                                       the second column')
    mutation_parser.add_argument('-gof', '--gene_order_file', required=True,
                                 help='Path to file containing tab-separated lists of genes on\
                                 each chromosme, in order of their position on the chromosome, one\
                                  chromosome per line')
    mutation_parser.add_argument('-b', '--bmr', type=float, required=True,
                                 help='Default background mutation rate')
    mutation_parser.add_argument('-bf', '--bmr_file',
                                 help='File listing gene-specific BMRs. If none, the default BMR\
                                       will be used for all genes.')
    return parser
    
def run(args):
    #if l not specified, set default based on test statistic 
    if not args.sizes:
        args.sizes = [5,10,15,20] if args.test_statistic == MAX_CC_SIZE else [3]
    
    #disallow finding delta by # of CCs of size >= l for HotNet2, since this is not currently
    #implemented correctly (and is non-trivial to implement)
    if not args.classic and args.test_statistic != MAX_CC_SIZE:
        raise ValueError("For HotNet2, the largest CC size test statistic must be used.")
    
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)

    if args.perm_type == "heat":
        infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
        addtl_genes = hnio.load_genes(args.permutation_genes_file) if args.permutation_genes_file else None
        deltas = get_deltas_for_heat(infmat, infmat_index, heat, addtl_genes, args.num_permutations,
                                     args.test_statistic, args.sizes, args.classic, args.parallel)
    elif args.perm_type == "mutations":
        infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
        deltas = get_deltas_for_mutations(args, infmat, infmat_index, heat_params)
    elif args.perm_type == "network":
        deltas = get_deltas_for_network(args.permuted_networks_path, heat, args.infmat_name,
                                         infmat_index, args.test_statistic, args.sizes,
                                         args.classic, args.num_permutations, args.parallel)
    else:
        raise ValueError("Invalid mutation permutation type: %s" % args.perm_type)
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump({"parameters": vars(args), "heat_parameters": heat_params,
               "deltas": deltas}, output_file, indent=4)
    if (args.output_file): output_file.close()

def get_deltas_for_network(permuted_networks_path, heat, infmat_name, index2gene, test_statistic,
                            sizes, classic, num_permutations, parallel):
    print "* Performing permuted network delta selection..."
    
    #construct list of paths to the first num_permutations     
    permuted_network_paths = [permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                              for i in range(1, num_permutations+1)]

    delta_selection_fn = delta.find_best_delta_by_largest_cc if test_statistic == "max_cc_size" \
                            else delta.find_best_delta_by_num_ccs 

    return delta.network_delta_selection(permuted_network_paths, infmat_name, index2gene, heat,
                                         sizes, not classic, parallel, delta_selection_fn)

def get_deltas_for_heat(infmat, index2gene, gene2heat, addtl_genes, num_permutations, test_statistic,
                        sizes, classic, parallel):
    print "* Performing permuted heat delta selection..."
    heat_permutations = permutations.permute_heat(gene2heat, num_permutations, addtl_genes, parallel)
    return get_deltas_from_heat_permutations(infmat, index2gene, heat_permutations, test_statistic,
                                             sizes, classic, parallel)

def get_deltas_for_mutations(args, infmat, index2gene, heat_params):
    print "* Performing permuted mutation data delta selection..."
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
    index2gene = hnio.load_index(args.infmat_index_file)
    
    heat_permutations = permutations.generate_mutation_permutation_heat(
                            heat_params["heat_fn"], heat_params["sample_file"],
                            heat_params["gene_file"], index2gene.values(), heat_params["snv_file"],
                            args.gene_length_file, args.bmr, args.bmr_file, heat_params["cna_file"],
                            args.gene_order_file, heat_params["cna_filter_threshold"],
                            heat_params["min_freq"], args.num_permutations, args.parallel)
    return get_deltas_from_heat_permutations(infmat, index2gene, heat_permutations, args.test_statistic,
                                             args.sizes, args.classic, args.parallel)

def get_deltas_from_heat_permutations(infmat, gene_index, heat_permutations, test_statistic, sizes,
                                      classic, parallel):
    delta_selection_fn = delta.find_best_delta_by_largest_cc if test_statistic == "max_cc_size" \
                            else delta.find_best_delta_by_num_ccs 
    
    deltas = delta.heat_delta_selection(infmat, gene_index, heat_permutations, sizes, not classic,
                                        parallel, delta_selection_fn)
    return deltas

if __name__ == "__main__": 
    run(get_parser().parse_args(sys.argv[1:]))
