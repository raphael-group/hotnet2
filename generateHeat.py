import sys
import json
import argparse
from hotnet import hnap, hnio, heat as hnheat

def parse_args(raw_args): 
    description = "Generates a JSON heat file for input to runHotNet2."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-o', '--output_file',
                               help='Output file.  If none given, output will be written to stdout.')
    
    subparsers = parser.add_subparsers(title='Heat score type')
    
    heat_parser = subparsers.add_parser('scores', help='Pre-computed heat scores', parents=[parent_parser])
    heat_parser.add_argument('-hf', '--heat_file', required=True,
                             help='Path to a tab-separated file containing a gene name in the first\
                                   column and the heat score for that gene in the second column of\
                                   each line.')
    heat_parser.add_argument('-ms', '--min_heat_score', type=float,
                             help='Minimum heat score for including genes in the resulting output\
                                   file. By default, all genes with positive heat scores will be\
                                   included.')
    heat_parser.add_argument('-gff', '--gene_filter_file', default=None,
                             help='Path to file listing genes whose heat scores should be\
                                   preserved, one per line. If present, all other heat scores\
                                   will be discarded.')
    heat_parser.add_argument('-e', '--excluded_genes_output_file',
                             help='File path to which the list of genes that were excluded from\
                                   the heat score output due to the specified filtering parameters\
                                   should be written, one gene per line. If no genes were filtered\
                                   and a path is specified, the resulting file will be empty.')
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
    
    return parser.parse_args(raw_args)

def valid_cna_filter_thresh(string):
        value = float(string)
        if value <= .5:
            raise argparse.ArgumentTypeError("cna_filter_threshold must be > .5")
        return value

def load_direct_heat(args):
    heat = hnio.load_heat_tsv(args.heat_file)
    
    #ensure that all heat scores are positive
    bad_genes = [gene for gene in heat if heat[gene] < 0]
    if bad_genes:
        raise ValueError("ERROR: All gene heat scores must be non-negative. There are %s genes with\
                          negative heat scores: %s" % (len(bad_genes), bad_genes))
    
    heat, excluded_genes, args.min_heat_score = hnheat.filter_heat(heat, args.min_heat_score)
    return heat, excluded_genes

def load_mutation_heat(args):
    samples = hnio.load_samples(args.sample_file) if args.sample_file else None
    genes = hnio.load_genes(args.gene_file) if args.gene_file else None
    snvs = hnio.load_snvs(args.snv_file, genes, samples)
    cnas = hnio.load_cnas(args.cna_file, genes, samples) if args.cna_file else []
    if args.cna_filter_threshold:
        cnas = hnheat.filter_cnas(cnas, args.cna_filter_threshold)
    
    if not samples:
        samples = set([snv.sample for snv in snvs] + [cna.sample for cna in cnas])
    return hnheat.mut_heat(len(samples), snvs, cnas, args.min_freq), None

def run(args):
    heat, heat_excluded_genes = args.heat_fn(args)
    
    filter_excluded_genes = []
    if args.heat_fn != load_mutation_heat and args.gene_filter_file:
        heat, filter_excluded_genes = hnheat.expr_filter_heat(heat,
                                                              hnio.load_genes(args.gene_filter_file))
    
    args.heat_fn = args.heat_fn.__name__
    output_dict = {"parameters": vars(args), "heat": heat}
    
    if args.heat_fn == "load_direct_heat":
        output_dict["excluded_genes"] = list(set().union(heat_excluded_genes, filter_excluded_genes))
        if args.excluded_genes_output_file:
            hnio.write_gene_list(args.excluded_genes_output_file, heat_excluded_genes)    
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump(output_dict, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
