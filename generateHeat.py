import argparse
import json
import sys
from hotnet2 import heat as hnheat, hnap, hnio

def get_parser():
    description = "Generates a JSON heat file for input to runHotNet2."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-o', '--output_file',
                               help='Output file.  If none given, output will be written to stdout.')
    parent_parser.add_argument('-n', '--name',
                               help='Name/Label describing the heat scores.')

    subparsers = parser.add_subparsers(title='Heat score type')

    heat_parser = subparsers.add_parser('scores', help='Pre-computed heat scores', parents=[parent_parser])
    heat_parser.add_argument('-hf', '--heat_file', required=True,
                             help='Path to a tab-separated file containing a gene name in the first\
                                   column and the heat score for that gene in the second column of\
                                   each line.')
    heat_parser.add_argument('-ms', '--min_heat_score', type=float, default=0,
                             help='Minimum heat score for genes to have their original heat score\
                                   in the resulting output file. Genes with score below this value\
                                   will be assigned score 0.')
    heat_parser.add_argument('-gff', '--gene_filter_file', default=None,
                             help='Path to file listing genes whose heat scores should be\
                                   preserved, one per line. If present, all other heat scores\
                                   will be discarded.')
    heat_parser.set_defaults(heat_fn=load_direct_heat)

    mutation_parser = subparsers.add_parser('mutation', help='Mutation data', parents=[parent_parser])
    mutation_parser.add_argument('--snv_file', required=True,
                                 help='Path to a tab-separated file containing SNVs where the first\
                                       column of each line is a sample ID and subsequent columns\
                                       contain the names of genes with SNVs in that sample. Lines\
                                       starting with "#" will be ignored.')
    mutation_parser.add_argument('--cna_file',
                                 help='Path to a tab-separated file containing CNAs where the first\
                                       column of each line is a sample ID and subsequent columns\
                                       contain gene names followed by "(A)" or "(D)" indicating an\
                                       amplification or deletion in that gene for the sample.\
                                       Lines starting with "#" will be ignored.')
    mutation_parser.add_argument('--sample_file', default=None,
                                 help='File listing samples. Any SNVs or CNAs in samples not listed\
                                       in this file will be ignored. If HotNet is run with mutation\
                                       permutation testing, all samples in this file will be eligible\
                                       for random mutations even if the sample did not have any\
                                       mutations in the real data. If not provided, the set of samples\
                                       is assumed to be all samples that are provided in the SNV\
                                       or CNA data.')
    mutation_parser.add_argument('--sample_type_file', default=None,
                                 help='File listing type (e.g. cancer, datasets, etc.) of samples\
                                       (see --sample_file). Each line is a space-separated row\
                                       listing one sample and its type. The sample types are used\
                                       for creating the HotNet(2) web output.')
    mutation_parser.add_argument('--gene_file', default=None,
                                 help='File listing tested genes. SNVs or CNAs in genes not listed\
                                       in this file will be ignored. If HotNet is run with mutation\
                                       permutation testing, every gene in this file will be eligible\
                                       for random mutations even if the gene did not have mutations\
                                       in any samples in the original data. If not provided, the set\
                                       of tested genes is assumed to be all genes that have mutations\
                                       in either the SNV or CNA data.')
    mutation_parser.add_argument('--min_freq', type=int, default=1,
                                 help='The minimum number of samples in which a gene must have an\
                                       SNV to be considered mutated in the heat score calculation.')
    mutation_parser.add_argument('--cna_filter_threshold', type=valid_cna_filter_thresh,
                                 default=None,
                                 help='Proportion of CNAs in a gene across samples that must share\
                                       the same CNA type in order for the CNAs to be included. This\
                                       must either be > .5, or the default, None, in which case all\
                                       CNAs will be included.')
    mutation_parser.set_defaults(heat_fn=load_mutation_heat)

    oncodrive_parser = subparsers.add_parser('oncodrive', help='Oncodrive scores', parents=[parent_parser])
    oncodrive_parser.add_argument('--fm_scores', required=True, help='Oncodrive-FM scores (gene to q-value).')
    oncodrive_parser.add_argument('--cis_amp_scores', required=True,
                                  help='Oncodrive-CIS scores (gene to q-value); amplifications only.')
    oncodrive_parser.add_argument('--cis_del_scores', required=True,
                                  help='Oncodrive-CIS scores (gene to q-value); deletions only.')
    oncodrive_parser.add_argument('--fm_threshold', type=float, default=0.2,
                                  help='Maximum Oncodrive-FM q-value threshold')
    oncodrive_parser.add_argument('--cis_threshold', type=float, default=0.2,
                                  help='Maximum Oncodrive-CIS q-value threshold')
    oncodrive_parser.add_argument('--cis', default=False, action='store_true',
                                  help='Flag whether to include Oncodrive-CIS scores when generating '\
                                        'the Oncodrive heat file.')
    oncodrive_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    oncodrive_parser.set_defaults(heat_fn=load_oncodrive_heat)

    mutsig_parser = subparsers.add_parser('mutsig', help='MutSig scores', parents=[parent_parser])
    mutsig_parser.add_argument('--mutsig_score_file', required=True, help='MutSig score file (gene to q-value).')
    mutsig_parser.add_argument('--threshold', type=float, default=1.0, help='Maximum q-value threshold.')
    mutsig_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    mutsig_parser.set_defaults(heat_fn=load_mutsig_heat)

    music_parser = subparsers.add_parser('music', help='MuSiC scores', parents=[parent_parser])
    music_parser.add_argument('--music_score_file', required=True, help='MuSiC score file (gene to q-value).')
    music_parser.add_argument('--threshold', type=float, default=1.0, help='Maximum q-value threshold.')
    music_parser.add_argument('--max_heat', type=float, default=15, help='Max heat')
    music_parser.add_argument('--gene_filter_file', default=None,
                                 help='File listing genes whose heat scores should be preserved.\
                                       If present, all other heat scores will be discarded.')
    music_parser.set_defaults(heat_fn=load_music_heat)

    return parser

def valid_cna_filter_thresh(string):
        value = float(string)
        if value <= .5:
            raise argparse.ArgumentTypeError("cna_filter_threshold must be > .5")
        return value

def load_direct_heat(args):
    heat = hnio.load_heat_tsv(args.heat_file)
    print "* Loading heat scores for %s genes" % len(heat)

    #ensure that all heat scores are positive
    bad_genes = [gene for gene in heat if heat[gene] < 0]
    if bad_genes:
        raise ValueError("ERROR: All gene heat scores must be non-negative. There are %s genes with\
                          negative heat scores: %s" % (len(bad_genes), bad_genes))

    heat, _ = hnheat.filter_heat(heat, args.min_heat_score, True,
                                 'Assigning score 0 to ## genes with score below %s' % args.min_heat_score)
    return heat

def load_mutation_heat(args):
    genes = hnio.load_genes(args.gene_file) if args.gene_file else None
    samples = hnio.load_samples(args.sample_file) if args.sample_file else None
    snvs = hnio.load_snvs(args.snv_file, genes, samples)
    cnas = hnio.load_cnas(args.cna_file, genes, samples) if args.cna_file else []
    if args.cna_filter_threshold:
        cnas = hnheat.filter_cnas(cnas, args.cna_filter_threshold)

    if not samples:
        samples = set([snv.sample for snv in snvs] + [cna.sample for cna in cnas])
    if not genes:
        genes = set([snv.gene for snv in snvs] + [cna.gene for cna in cnas])
    return hnheat.mut_heat(genes, len(samples), snvs, cnas, args.min_freq)

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
    if args.heat_fn != load_mutation_heat and args.gene_filter_file:
        heat = hnheat.reconcile_heat_with_tested_genes(heat, hnio.load_genes(args.gene_filter_file))

    args.heat_fn = args.heat_fn.__name__
    output_dict = {"parameters": vars(args), "heat": heat}

    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
