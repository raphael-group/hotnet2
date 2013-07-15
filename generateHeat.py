import hnap
import hnio
import heat
import sys
import json

def parse_args(raw_args): 
    description = "Generates a JSON heat file for input to runHotNet2."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    
    subparsers = parser.add_subparsers(title='Heat score type')
    
    heat_parser = subparsers.add_parser('direct', help='Direct heat scores')
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
    return hnio.load_heat_tsv(args.heat_file)

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
    heat = args.heat_fn(args)
    if args.gene_filter_file:
        heat = heat.expr_filter_heat(heat, hnio.load_gene_list(args.gene_filter_file))
    
    args.heat_fn = args.heat_fn.__name__
    output_dict = {"parameters": vars(args), "heat": heat}
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
