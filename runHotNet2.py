# -*- coding: iso-8859-1 -*-
import sys
import json
import os
import scipy.io
from hotnet import hnap, hnio, hotnet2 as hn, permutations, stats
from hotnet.constants import *

def parse_args(raw_args): 
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
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=3,
                        help='Minimum size connected components that should be returned.')
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
    parent_parser.add_argument('--parallel', dest='parallel', action='store_true',
                               help='Run permutation tests in parallel. Only recommended for machines\
                                     with at least 8 cores.')
    parent_parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                               help='Run permutation tests sequentially. Recommended for machines\
                                     with fewer than 8 cores.')
    parent_parser.set_defaults(parallel=False)
    
    subparsers = parser.add_subparsers(title='Heat score type', dest='permutation_type')
    
    subparsers.add_parser('none', help='Do not perform statistical significance permutation tests')
    
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-pgf', '--permutation_genes_file',
                             help='Path to file containing a list of additional genes that can have\
                                   permuted heat values assigned to them in permutation tests')
    
    precomp_parser = subparsers.add_parser('precomputed', help='Use precomputed datasets',
                                           parents=[parent_parser])
    precomp_parser.add_argument('-dp', '--datasets_path', required=True,
                                help='Path to datasets to use for significance testing. Include ' +
                                      ITERATION_REPLACEMENT_TOKEN + ' in the path to be replaced\
                                      with the iteration number.')
    
    return parser.parse_args(raw_args)

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    if len(os.listdir(args.output_directory)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")
    
    # load data
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)
  
    # compute similarity matrix and extract connected components
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h)
    
    # only calculate permuted data sets for significance testing once
    if args.permutation_type != "none":
        if args.permutation_type == "heat":
            print "* Generating heat permutations for statistical significance testing" 
            extra_genes = hnio.load_genes(args.permutation_genes_file) if args.permutation_genes_file \
                            else None
            heat_permutations = permutations.permute_heat(heat, args.num_permutations, extra_genes,
                                                          args.parallel)
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
        
        G = hn.weighted_graph(sim, gene_index, delta)
        ccs = hn.connected_components(G, args.min_cc_size)
        
        # calculate significance
        if args.permutation_type != "none":
            sizes2stats = calculate_significance(args, infmat, infmat_index, G, delta, heat_permutations)
        
        #sort ccs list such that genes within components are sorted alphanumerically, and components
        #are sorted first by length, then alphanumerically by name of the first gene in the component 
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)
        
        #write output
        hnio.write_components_as_tsv(os.path.abspath(delta_out_dir) + "/" + COMPONENTS_TSV, ccs)
        args.delta = delta
        output_dict = {"parameters": vars(args), "heat_parameters": heat_params,
                       "sizes": hn.component_sizes(ccs), "components": ccs}
        if args.permutation_type != "none":
            output_dict["statistics"] = sizes2stats
            hnio.write_significance_as_tsv(os.path.abspath(delta_out_dir) + "/" + SIGNIFICANCE_TSV,
                                           sizes2stats)
        
        json_out = open(os.path.abspath(delta_out_dir) + "/" + JSON_OUTPUT, 'w')
        json.dump(output_dict, json_out, indent=4)
        json_out.close()

def calculate_significance(args, infmat, infmat_index, G, delta, heat_permutations):
    sizes = range(args.cc_start_size, args.cc_stop_size+1)
    
    print "\t- Using no. of components >= k (k \\in",
    print "[%s, %s]) as statistic" % (min(sizes), max(sizes))

    #size2counts is dict(size -> (list of counts, 1 per permutation))
    sizes2counts = stats.calculate_permuted_cc_counts(infmat, infmat_index, heat_permutations,
                                                      delta, sizes, args.parallel)
    real_counts = stats.num_components_min_size(G, sizes)
    size2real_counts = dict(zip(sizes, real_counts))
    return stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
