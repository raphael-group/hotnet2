# -*- coding: iso-8859-1 -*-
import hnap
import hnio
import hotnet2 as hn
import permutations
import stats
import sys
import scipy.io
import json

def parse_args(raw_args):  
    
    description = "Runs generalized HotNet2.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python runHotnet2.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-mn', '--infmat_name', required=True, default='Li',
                        help='Variable name of the influence matrix in the .mat file')
    parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                        help='Gene-index file for the influence matrix.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='Heat score file')
    parser.add_argument('-d', '--delta', type=float, required=True,
                        help='Weight threshold for edge removal')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=3,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-c', '--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    parser.add_argument('-n', '--num_permutations', type=int, required=True,
                        help='Number of permutation tests to run; set to 0 to skip running\
                              permutation tests.')
    parser.add_argument('-pgf', '--permutation_genes_file',
                        help='Path to file containing a list of additional genes that can have\
                              permuted heat values assigned to them in permutation tests') 
    parser.add_argument('-s', '--cc_start_size', type=int, default=2,
                        help='Smallest connected component size to count')
    parser.add_argument('-l', '--cc_stop_size', type=int, default=10,
                        help='Largest connected component size to count')
    parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    parser.add_argument('-p', '--parallel', default=False, action='store_true',
                        help='Include flag to run permutation tests in parallel.')
    
    return parser.parse_args(raw_args)

def run(args):
    # load data
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
    infmat_index = hnio.load_index(args.infmat_index_file)  #dict from indices to gene names 
    heat = hnio.load_heat(args.heat_file)                   #dict from gene names to heat scores
  
    # compute similarity matrix and extract connected components
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h, not args.classic)
    G = hn.weighted_graph(sim, gene_index, args.delta, not args.classic)
    ccs = hn.connected_components(G, args.min_cc_size)
    
    # calculate significance
    if args.num_permutations > 0:
        extra_genes = hnio.load_gene_list(args.permutation_genes_file)
        heat_permutations = permutations.permute_heat(heat, args.num_permutations, extra_genes,
                                                      args.parallel)
        sizes = range(args.cc_start_size, args.cc_stop_size+1)
    
        #size2counts is dict(size -> (list of counts, 1 per permutation))
        sizes2counts = stats.calculate_permuted_cc_counts(infmat, infmat_index, heat_permutations,
                                                          args.delta, sizes, not args.classic,
                                                          args.parallel)
        real_counts = stats.num_components_min_size(G, sizes)
        size2real_counts = dict(zip(sizes, real_counts))
        sizes2stats = stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)
    
    # write output
    output_dict = {"parameters": vars(args), "sizes": hn.component_sizes(ccs), "components": ccs}
    if args.num_permutations > 0:
        output_dict["statistics"] = sizes2stats
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run( parse_args(sys.argv[1:]) )
