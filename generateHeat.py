import hnap
import hnio
import heat as hnheat
import sys
import json
import argparse

def parse_args(raw_args): 
    description = "Generates a JSON heat file for input to runHotNet2."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-o', '--output_file',
                               help='Output file.  If none given, output will be written to stdout.')
    
    subparsers = parser.add_subparsers(title='Heat score type')
    
    heat_parser = subparsers.add_parser('direct', help='Direct heat scores', parents=[parent_parser])
    heat_parser.add_argument('-hf', '--heat_file', required=True, help='Heat score file')
    heat_parser.add_argument('--gene_filter_file', default=None,
                             help='File listing genes whose heat scores should be preserved.\
                                   If present, all other heat scores will be discarded.')
    heat_parser.set_defaults(heat_fn=load_direct_heat)
    
    mutation_parser = subparsers.add_parser('mutation', help='Mutation data', parents=[parent_parser])
    mutation_parser.add_argument('--snv_file', required=True, help='SNV file')
    mutation_parser.add_argument('--cna_file', required=True, help='CNA file')
    mutation_parser.add_argument('--sample_file',
                                 help='File listing samples. Any SNVs or CNAs in samples not listed\
                                       in this file will be ignored. If HotNet is run with mutation\
                                       permutation testing, all samples in this file will be eligible\
                                       for random mutations even if the sample did not have any\
                                       mutations in the real data. If not provided, the set of samples\
                                       is assumed to be all samples that are provided in the SNV\
                                       or CNA data.')
    mutation_parser.add_argument('--gene_file',
                                 help='File listing tested genes. SNVs or CNAs in genes not listed\
                                       in this file will be ignored. If HotNet is run with mutation\
                                       permutation testing, every gene in this file will be eligible\
                                       for random mutations even if the gene did not have mutations\
                                       in any samples in the original data. If not provided, the set\
                                       of tested genes is assumed to be all genes that have mutations\
                                       in either the SNV or CNA data.')
    mutation_parser.add_argument('--min_freq', type=int, default=1, help='Minimum frequency')
    mutation_parser.add_argument('--cna_filter_threshold', type=valid_cna_filter_thresh,
                                 default=None,
                                 help='Proportion of CNAs in a gene across samples that must share\
                                       the same CNA type in order for the CNAs to be included. This\
                                       must either be > .5, or the default, None, in which case all\
                                       CNAs will be included')
    mutation_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    mutation_parser.set_defaults(heat_fn=load_mutation_heat)
    
    oncodrive_parser = subparsers.add_parser('oncodrive', help='Oncodrive scores', parents=[parent_parser])
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
    
    mutsig_parser = subparsers.add_parser('mutsig', help='MutSig scores', parents=[parent_parser])
    mutsig_parser.add_argument('--mutsig_score_file', required=True, help='MutSig score file')
    mutsig_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    mutsig_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    mutsig_parser.set_defaults(heat_fn=load_mutsig_heat)
    
    music_parser = subparsers.add_parser('music', help='MuSiC scores', parents=[parent_parser])
    music_parser.add_argument('--music_score_file', required=True, help='MuSiC score file')
    music_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    music_parser.add_argument('--max_heat', type=float, default=15, help='Max heat')
    music_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    music_parser.set_defaults(heat_fn=load_music_heat)
    
    return parser.parse_args(raw_args)

def valid_cna_filter_thresh(string):
        value = float(string)
        if value <= .5:
            raise argparse.ArgumentTypeError("cna_filter_threshold must be > .5")
        return value

def load_direct_heat(args):
    return hnio.load_heat_tsv(args.heat_file)

def load_mutation_heat(args):
    samples = hnio.load_samples(args.sample_file) if args.sample_file else None
    genes = hnio.load_genes(args.gene_file) if args.gene_file else None
    snvs = hnio.load_snvs(args.snv_file, genes, samples)
    cnas = hnio.load_cnas(args.cna_file, genes, samples)
    if args.cna_filter_threshold:
        cnas = hnheat.filter_cnas(cnas, args.cna_filter_threshold)
    
    if not samples:
        samples = set([snv.sample for snv in snvs] + [cna.sample for cna in cnas])
    return hnheat.mut_heat(len(samples), snvs, cnas, args.min_freq)

def load_oncodrive_heat(args):
    gene2heat = hnio.load_oncodrive_data(args.fm_scores, args.cis_amp_scores, args.cis_del_scores)
    return hnheat.fm_heat(gene2heat, args.fm_threshold, args.cis_threshold, args.cis)
    
def load_mutsig_heat(args):
    gene2mutsig = hnio.load_mutsig_scores(args.mutsig_score_file)
    return hnheat.mutsig_heat(gene2mutsig, args.threshold)

def load_music_heat(args):
    gene2music = hnio.load_music_scores(args.music_score_file)
    return hnheat.music_heat(gene2music, args.threshold, args.max_heat)

def run(args):
    heat = args.heat_fn(args)
    if args.gene_filter_file:
        heat = hnheat.expr_filter_heat(heat, hnio.load_genes(args.gene_filter_file))
    
    args.heat_fn = args.heat_fn.__name__
    output_dict = {"parameters": vars(args), "heat": heat}
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
