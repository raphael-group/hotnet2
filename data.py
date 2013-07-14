import hnio
import hnap
import sys
from math import log10
import scipy

def parse_args(raw_args):
    description = "Data module"
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    subparsers = parser.add_subparsers(title='Heat scores')

    #create subparser for options for permuting networks
    mutation_parser = subparsers.add_parser('mutation', help='Mutation data')
    mutation_parser.add_argument('--snv_file', required=True, help='SNV file')
    mutation_parser.add_argument('--cna_file', required=True, help='CNA file')
    mutation_parser.add_argument('--sample_file', required=True, help='Sample file')
    mutation_parser.add_argument('--gene_file', required=True, help='Gene file')
    mutation_parser.add_argument('--min_freq', type=int, default=11, help="Minimum frequency")
    mutation_parser.add_argument('--output_file', required=True, help='Output file')
    mutation_parser.set_defaults(heat_fn=mutation)
    
    
#ONCODRIVE_DIR = "/data/compbio/datasets/TCGA/PanCancer/UPF-Barcelona/oncodrive"
#oncodrive_fm  = "/data/compbio/datasets/TCGA/PanCancer/UPF-Barcelona/oncodrive/pan12_oncoFM_z_values_abr2013_zscores.tsv" % ONCODRIVE_DIR
#oncodrive_amp = "/data/compbio/datasets/TCGA/PanCancer/UPF-Barcelona/oncodrive/cis/2013-06-06/filtered/pan12_OncodriveCIS_Over.txt" % ONCODRIVE_DIR 
#oncodrive_del = "/data/compbio/datasets/TCGA/PanCancer/UPF-Barcelona/oncodrive/cis/2013-06-06/filtered/pan12_OncodriveCIS_Under.txt" % ONCODRIVE_DIR
#GENE2ONCO = load_oncodrive_data(oncodrive_fm, oncodrive_amp, oncodrive_del)
#genes, heat = fm_heat(GENE2ONCO, 0.2, 0.2, CIS=True)
#     gene2heat = hnio.load_oncodrive_data(args.fm_scores, args.cis_amp_scores, args.cis_del_scores)
#     return fm_heat(gene2heat, args.fm_threshold, args.cis_threshold, args.CIS)
    
    oncodrive_parser = subparsers.add_parser('oncodrive', help='Oncodrive scores')
    oncodrive_parser.add_argument('--fm_scores', required=True, help='???')
    oncodrive_parser.add_argument('--cis_amp_scores', required=True, help='???')
    oncodrive_parser.add_argument('--cis_del_scores', required=True, help='???')
    oncodrive_parser.add_argument('--fm_threshold', type=float, default=0.2, help='???')
    oncodrive_parser.add_argument('--cis_threshold', type=float, default=0.2, help='???')
    oncodrive_parser.add_argument('--cis', default=False, action='store_true', help='???')
    oncodrive_parser.add_argument('--output_file', required=True, help='Output file')
    oncodrive_parser.set_defaults(heat_fn=oncodrive)
    
    mutsig_parser = subparsers.add_parser('mutsig', help='MutSig scores')
    mutsig_parser.add_argument('--mutsig_score_file', required=True, help='MutSig score file')
    mutsig_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    mutsig_parser.add_argument('--output_file', required=True, help='Output file')
    mutsig_parser.set_defaults(heat_fn=mutsig)
    
    music_parser = subparsers.add_parser('music', help='MuSiC scores')
    music_parser.add_argument('--music_score_file', required=True, help='MuSiC score file')
    music_parser.add_argument('--threshold', type=float, default=1.0, help='Threshold...no idea what this is')
    music_parser.add_argument('--max_heat', type=float, default=15, help='Max heat')
    music_parser.add_argument('--output_file', required=True, help='Output file')
    music_parser.set_defaults(heat_fn=music)
    
    return parser.parse_args(raw_args)

def num_snvs(mutation_list):
    return len([ p for p, muts in mutation_list.items() if "snv" in muts ])

def num_cnas(mutation_list):
    return len([ p for p, muts in mutation_list.items() if "amp" in muts or "del" in muts ])

def mut_heat((genes, samples, gene2mutations), min_freq):
    print "* Creating mutation fraction heat map..."
    n = float(len(samples))
    heat = dict([(g, len( heat ) / n) for g, heat in gene2mutations.items()
                 if num_snvs(heat) > min_freq or num_cnas(heat) > 0])
    genes = sorted( heat.keys() )
    samples = set([ s for g in genes for s in gene2mutations[g] ])
    print "\t- Including", len(genes), "genes at min frequency", min_freq
    print "\t- Samples with no mutations (not taken into account for computing",
    print "heat):", int(n) - len(samples)
    return genes, heat

NULL = 100
def fm_heat(gene2heat, fm_threshold, cis_threshold=0.01, CIS=False):
    print "* Creating oncodrive heat map..."
    if CIS: print "\tIncluding CIS scores at threshold", cis_threshold, "..."
    heat = dict()
    src_fm, src_cis_amp, src_cis_del = 0, 0, 0
    for g, scores in gene2heat.items():
        if CIS:
            del_score = scores["del"] if scores["del"] < cis_threshold else NULL
            amp_score = scores["amp"] if scores["amp"] < cis_threshold else NULL
            fm_score  = scores["fm"] if scores["fm"] < fm_threshold else NULL
            if fm_score == NULL and amp_score == NULL and del_score == NULL: continue
            min_val = min(del_score, amp_score, fm_score)
            heat[g] = -log10( min_val )
            if min_val == scores["fm"]: src_fm += 1
            elif min_val == scores["amp"]: src_cis_amp += 1
            elif min_val == scores["del"]: src_cis_del += 1
        else:
            if scores["fm"] >= fm_threshold: continue
            heat[g] = -log10(scores["fm"])
            src_fm += 1
    print "\t- Genes using FM score:", src_fm
    print "\t- Genes using CIS AMP score:", src_cis_amp
    print "\t- Genes using CIS DEL score:", src_cis_del
 
    return sorted(heat.keys()), heat

def mutsig_heat(gene2mutsig, threshold=1.0):
    print "* Creating MutSig heat map..."
    gene2heat = dict([(gene, -log10(score["qval"]))
                      for gene, score in gene2mutsig.items()
                      if score["qval"] < threshold])
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat.keys(), gene2heat

def music_heat(gene2music, threshold=1.0, max_heat=15):
    print "* Creating MuSiC heat map..."
    print "\t- Mapping q-values of 0 to", max_heat
    def music_heat(qvals):
        heat = scipy.median([ qvals["FDR_CT"], qvals["FDR_LRT"], qvals["FDR_FCPT"] ])
        return -log10(heat) if heat != 0 else max_heat
    gene2heat = dict([(gene, music_heat(scores)) for gene, scores in gene2music.items()
                      if scipy.median(scores.values()) < threshold])
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat.keys(), gene2heat
 
def mutation(args):
    genes_samples_gene2heat, samples_w_tys = hnio.load_mutation_data(args.snv_file, args.cna_file, args.sample_file, args.gene_file)
    return mut_heat(genes_samples_gene2heat, args.min_freq)

def oncodrive(args):
    gene2heat = hnio.load_oncodrive_data(args.fm_scores, args.cis_amp_scores, args.cis_del_scores)
    return fm_heat(gene2heat, args.fm_threshold, args.cis_threshold, args.cis)
    
def mutsig(args):
    gene2mutsig = hnio.load_mutsig_scores(args.mutsig_score_file)
    return mutsig_heat(gene2mutsig, args.threshold)

def music(args):
    gene2music = hnio.load_music_scores(args.music_score_file)
    return music_heat(gene2music, args.threshold, args.max_heat)

def run(args):
    genes, heat = args.heat_fn(args)
    hnio.save_heat(heat, args.output_file)
    
#  
# def expr_filter_heat( gene2heat ):
#     gene2heat = dict([(g, h) for g, h in gene2heat.items() if g in EXPRESSED_GENES])
#     return gene2heat.keys(), gene2heat

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
