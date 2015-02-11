# -*- coding: iso-8859-1 -*-
import json
import os
import sys
import os.path
sys.path.append(os.path.split(os.path.split(sys.argv[0])[0])[0])
import numpy as np
from hotnet2 import hnap, hnio, hotnet2 as hn, permutations as p, stats
from hotnet2.constants import ITERATION_REPLACEMENT_TOKEN, JSON_OUTPUT, COMPONENTS_TSV, SIGNIFICANCE_TSV

def get_parser():
    description = "Runs generalized HotNet2.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python runHotnet2.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-mn', '--infmat_name', default='PPR',
                        help='Variable name of the influence matrix in the .mat file')
    parser.add_argument('-if', '--infmat_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='JSON heat score file generated via generateHeat.py')
    parser.add_argument('-d', '--deltas', nargs='+', type=float, required=True,
                        help='Weight threshold for edge removal')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=2,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-c', '--classic', default=False, action='store_true',
                        help='Run classic HotNet (rather than HotNet2).')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory. Files results.json, components.txt, and\
                              significance.txt will be generated in subdirectories for each delta.')
    
    #parent parser for arguments common to all permutation types
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                               help='Number of permutation tests to run')
    parent_parser.add_argument('-s', '--cc_start_size', type=int, default=2,
                               help='Smallest connected component size to count in permutation tests')
    parent_parser.add_argument('-l', '--cc_stop_size', type=int, default=10,
                               help='Largest connected component size to count in permutation tests')
    parent_parser.add_argument('-c', '--num_cores', type=int, default=1,
                               help='Number of cores to use for running permutation tests in\
                               parallel. If -1, all available cores will be used.')
    
    subparsers = parser.add_subparsers(title='Permutation type', dest='permutation_type')
    
    subparsers.add_parser('none', help='Do not perform statistical significance permutation tests')
    
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-pgf', '--permutation_genes_file',
                             help='Path to file containing a list of additional genes that can have\
                                   permuted heat values assigned to them in permutation tests')
    
    mutation_parser = subparsers.add_parser('mutations', help='Permute mutation data',
                                             parents=[parent_parser])
    mutation_parser.add_argument('-glf', '--gene_length_file', required=True,
                                 help='Path to tab-separated file containing gene names in the\
                                       first column and the length of the gene in base pairs in\
                                       the second column')
    mutation_parser.add_argument('-gof', '--gene_order_file', required=True,
                                 help='Path to file containing tab-separated lists of genes on\
                                 each chromosome, in order of their position on the chromosome,\
                                 one chromosome per line')
    mutation_parser.add_argument('-b', '--bmr', type=float, required=True,
                                 help='Default background mutation rate')
    mutation_parser.add_argument('-bf', '--bmr_file',
                                 help='File listing gene-specific BMRs. If none, the default BMR\
                                       will be used for all genes.')
    
    #create subparser for options for permuting networks
    network_parser = subparsers.add_parser('network', help='Permute networks', parents=[parent_parser])
    network_parser.add_argument('-pnp', '--permuted_networks_path', required=True,
                                help='Path to influence matrices for permuted networks.\
                                      Include ' + ITERATION_REPLACEMENT_TOKEN + ' in the\
                                      path to be replaced with the iteration number')
    
    precomp_parser = subparsers.add_parser('precomputed', help='Use precomputed datasets',
                                           parents=[parent_parser])
    precomp_parser.add_argument('-dp', '--datasets_path', required=True,
                                help='Path to datasets to use for significance testing. Include ' +
                                      ITERATION_REPLACEMENT_TOKEN + ' in the path to be replaced\
                                      with the iteration number.')
    
    return parser

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    if len(os.listdir(args.output_directory)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")
    
    # load data
    infmat = hnio.load_hdf5(args.infmat_file)[args.infmat_name]
    full_index2gene = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)
  
    # compute similarity matrix
    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, not args.classic)
    
    # only calculate permuted data sets for significance testing once
    if args.permutation_type != "none":
        if args.permutation_type == "heat":
            print "* Generating heat permutations for statistical significance testing" 
            extra_genes = hnio.load_genes(args.permutation_genes_file) \
                            if args.permutation_genes_file else None
            heat_permutations = p.permute_heat(heat, full_index2gene.values(),
                                               args.num_permutations, extra_genes, args.num_cores)
        elif args.permutation_type == "mutations":
            if heat_params["heat_fn"] != "load_mutation_heat":
                    raise RuntimeError("Heat scores must be based on mutation data to perform\
                                        significance testing based on mutation data permutation.")
            print "* Generating mutation permutations for statistical significance testing"
            heat_permutations = p.generate_mutation_permutation_heat(
                                    heat_params["heat_fn"], heat_params["sample_file"],
                                    heat_params["gene_file"], full_index2gene.values(),
                                    heat_params["snv_file"], args.gene_length_file, args.bmr,
                                    args.bmr_file, heat_params["cna_file"], args.gene_order_file,
                                    heat_params["cna_filter_threshold"], heat_params["min_freq"],
                                    args.num_permutations, args.num_cores)
        elif args.permutation_type == "network":
            pass    #nothing to do right now
        elif args.permutation_type == "precomputed":
            heat_file_paths = [args.datasets_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                               for i in range(1, args.num_permutations+1)]
            heat_permutations = [hnio.load_heat_tsv(heat_file) for heat_file in heat_file_paths]
        else:
            raise ValueError("Unrecognized permutation type %s" % (args.permutation_type))
    
    for delta in args.deltas:
        delta_out_dir = args.output_directory + "/delta_" + str(delta)
        if not os.path.isdir(delta_out_dir):
            os.mkdir(delta_out_dir)
        
        G = hn.weighted_graph(sim, index2gene, delta, not args.classic)
        ccs = hn.connected_components(G, args.min_cc_size)
        
        # calculate significance
        if args.permutation_type != "none":
            if args.permutation_type == "network":
                sizes2stats = calculate_significance_network(args, args.permuted_networks_path,
                                                             full_index2gene, G, heat, delta,
                                                             args.num_permutations)
            else:
                sizes2stats = calculate_significance(args, infmat, full_index2gene, G, delta,
                                                     heat_permutations)
        
        #sort ccs list such that genes within components are sorted alphanumerically, and components
        #are sorted first by length, then alphanumerically by name of the first gene in the component 
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)
        
        #write output
        hnio.write_components_as_tsv(os.path.abspath(delta_out_dir) + "/" + COMPONENTS_TSV, ccs)
        args.delta = delta  # include delta in parameters section of output JSON
        output_dict = {"parameters": vars(args), "heat_parameters": heat_params,
                       "sizes": hn.component_sizes(ccs), "components": ccs}
        if args.permutation_type != "none":
            output_dict["statistics"] = sizes2stats
            hnio.write_significance_as_tsv(os.path.abspath(delta_out_dir) + "/" + SIGNIFICANCE_TSV,
                                           sizes2stats)
        
        json_out = open(os.path.abspath(delta_out_dir) + "/" + JSON_OUTPUT, 'w')
        json.dump(output_dict, json_out, indent=4)
        json_out.close()

def calculate_significance_network(args, permuted_networks_path, index2gene, G, heat, delta, num_permutations):
    sizes = range(args.cc_start_size, args.cc_stop_size+1)
    
    print "\t- Using no. of components >= k (k \\in",
    print "[%s, %s]) as statistic" % (min(sizes), max(sizes))
    
    permuted_network_paths = [permuted_networks_path.replace(ITERATION_REPLACEMENT_TOKEN, str(i))
                              for i in range(1, num_permutations+1)]

    #size2counts is dict(size -> (list of counts, 1 per permutation))
    sizes2counts = stats.calculate_permuted_cc_counts_network(permuted_network_paths, args.infmat_name,
                                                        index2gene, heat, delta, sizes,
                                                        not args.classic, args.num_cores)
    
    real_counts = stats.num_components_min_size(G, sizes)
    size2real_counts = dict(zip(sizes, real_counts))
    return stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)

def calculate_significance(args, infmat, infmat_index, G, delta, heat_permutations):
    sizes = range(args.cc_start_size, args.cc_stop_size+1)
    
    print "\t- Using no. of components >= k (k \\in",
    print "[%s, %s]) as statistic" % (min(sizes), max(sizes))

    #size2counts is dict(size -> (list of counts, 1 per permutation))
    sizes2counts = stats.calculate_permuted_cc_counts(infmat, infmat_index, heat_permutations,
                                                      delta, sizes, not args.classic,
                                                      args.num_cores)
    real_counts = stats.num_components_min_size(G, sizes)
    size2real_counts = dict(zip(sizes, real_counts))
    return stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)

if __name__ == "__main__": 
    run(get_parser().parse_args(sys.argv[1:]))
