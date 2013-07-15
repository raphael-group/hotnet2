# -*- coding: iso-8859-1 -*-
import hnap
import hnio
import heat
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
    
    subparsers = parser.add_subparsers(title='Heat scores')
    
    heat_parser = subparsers.add_parser('heat', help='Direct heat scores')
    heat_parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    heat_parser.add_argument('--gene_filter_file', default=None,
                             help='File listing genes whose heat scores should be preserved.\
                                   If present, all other heat scores will be discarded.')
    heat_parser.set_defaults(heat_fn=load_direct_heat)
    
    mutation_parser = subparsers.add_parser('mutation', help='Mutation data')
    mutation_parser.add_argument('--snv_file', required=True, help='SNV file')
    mutation_parser.add_argument('--cna_file', required=True, help='CNA file')
    mutation_parser.add_argument('--sample_file', required=True, help='Sample file')
    mutation_parser.add_argument('--gene_file', required=True, help='Gene file')
    mutation_parser.add_argument('--min_freq', type=int, default=11, help='Minimum frequency')
    mutation_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    mutation_parser.set_defaults(heat_fn=load_mutation_heat)
    
    oncodrive_parser = subparsers.add_parser('oncodrive', help='Oncodrive scores')
    oncodrive_parser.add_argument('--fm_scores', required=True, help='???')
    oncodrive_parser.add_argument('--cis_amp_scores', required=True, help='???')
    oncodrive_parser.add_argument('--cis_del_scores', required=True, help='???')
    oncodrive_parser.add_argument('--fm_threshold', type=float, default=0.2, help='???')
    oncodrive_parser.add_argument('--cis_threshold', type=float, default=0.2, help='???')
    oncodrive_parser.add_argument('--cis', default=False, action='store_true', help='???')
    oncodrive_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    oncodrive_parser.set_defaults(heat_fn=load_oncodrive_heat)
    
    mutsig_parser = subparsers.add_parser('mutsig', help='MutSig scores')
    mutsig_parser.add_argument('--mutsig_score_file', required=True, help='MutSig score file')
    mutsig_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    mutsig_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    mutsig_parser.set_defaults(heat_fn=load_mutsig_heat)
    
    music_parser = subparsers.add_parser('music', help='MuSiC scores')
    music_parser.add_argument('--music_score_file', required=True, help='MuSiC score file')
    music_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    music_parser.add_argument('--max_heat', type=float, default=15, help='Max heat')
    music_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    music_parser.set_defaults(heat_fn=load_music_heat)
    
    return parser.parse_args(raw_args)

def load_direct_heat(args):
    return hnio.load_heat(args.heat_file)

def load_mutation_heat(args):
    genes_samples_gene2heat, _ = hnio.load_mutation_data(args.snv_file, args.cna_file, args.sample_file, args.gene_file)
    return heat.mut_heat(genes_samples_gene2heat, args.min_freq)

def load_oncodrive_heat(args):
    gene2heat = hnio.load_oncodrive_data(args.fm_scores, args.cis_amp_scores, args.cis_del_scores)
    return heat.fm_heat(gene2heat, args.fm_threshold, args.cis_threshold, args.cis)
    
def load_mutsig_heat(args):
    gene2mutsig = hnio.load_mutsig_scores(args.mutsig_score_file)
    return heat.mutsig_heat(gene2mutsig, args.threshold)

def load_music_heat(args):
    gene2music = hnio.load_music_scores(args.music_score_file)
    return heat.music_heat(gene2music, args.threshold, args.max_heat)

def run(args):
    # load data
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
    infmat_index = hnio.load_index(args.infmat_index_file)  #dict from indices to gene names 
            
    heat = args.heat_fn(args)                               #dict from gene names to heat scores
    if args.gene_filter_file:
        heat = heat.expr_filter_heat(heat, hnio.load_gene_list(args.gene_filter_file))
  
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
    args.heat_fn = args.heat_fn.__name__
    output_dict = {"parameters": vars(args), "sizes": hn.component_sizes(ccs), "components": ccs}
    if args.num_permutations > 0:
        output_dict["statistics"] = sizes2stats
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run( parse_args(sys.argv[1:]) )
