import json
from collections import defaultdict
from constants import *

################################################################################
# Data loading functions

def load_index(index_file):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_ppi_edges(edge_list_file):
    """Load PPI edges from file and return as a set of 2-tuples of gene indices.
    
    Arguments:
    edge_list_file -- path to file containing edges with a gene name in each of the first two
                      columns.
                      
    Note that edges are undirected, but each edge is represented as a single tuple in the
    returned set. Thus, to check whether a given pair of proteins interact, one must check
    for the presence of either ordered tuple.
    """
    arrs = [l.split() for l in open(edge_list_file)]
    return set([(int(arr[0]), int(arr[1])) for arr in arrs])

def load_heat_json(heat_file):
    with open(heat_file) as f:
        blob = json.load(f)
        return blob["heat"], blob["parameters"]

def load_heat_tsv(heat_file):
    arrs = [l.split() for l in open(heat_file)]
    return dict([(arr[0], float(arr[1])) for arr in arrs])

def load_genes(gene_file):
    """Load tested genes from a file and return as a set.
    
    Arguments:
    gene_file -- path to file containing gene names, one per line.
    
    """
    return set([l.strip() for l in open(gene_file)])

def load_gene_lengths(gene_lengths_file):
    arrs = [l.split() for l in open(gene_lengths_file)]
    return dict([(arr[0], int(arr[1])) for arr in arrs])

def load_gene_order(gene_order_file):
    """Load gene order file and return gene->chromosome and chromosome->ordered gene list mappings.
    
    Arguments:
    gene_order_file -- path to file containing tab-separated lists of genes on each chromosme,
                       one chromosome per line
    
    Note that numeric chromosome identifier used is simply the line number for the chromosome in
    the given file and does not indicate the true chromosome number.  
    """
    chromo2genes = {}
    gene2chromo = {}
    
    cid = 0
    for line in open(gene_order_file):
        genes = line.split()
        chromo2genes[cid] = genes
        gene2chromo.update([(gene, cid) for gene in genes])
        cid += 1
        
    return gene2chromo, chromo2genes

def load_gene_specific_bmrs(bmr_file):
    arrs = [l.split() for l in open(bmr_file)]
    return dict([(arr[0], float(arr[1])) for arr in arrs])  

def load_samples(sample_file):
    """Load sample IDs from a file and return as a set.
    
    Arguments:
    sample_file -- path to TSV file containing sample IDs as the first column. Any other columns
                   will be ignored.
    
    """
    return set([l.rstrip().split()[0] for l in open(sample_file)])

def include(item, whitelist):
    return item in whitelist if whitelist else True

#Returns list of Mutation objects representing SNVs
def load_snvs(snv_file, gene_wlst=None, sample_wlst=None):
    arrs = [l.rstrip().split("\t") for l in open(snv_file) if not l.startswith("#")]
    return [Mutation(arr[0], gene, SNV) for arr in arrs if include(arr[0], sample_wlst)
            for gene in arr[1:] if include(gene, gene_wlst)]

#Returns list of Mutation objects representing CNAs
def load_cnas(cna_file, gene_wlst=None, sample_wlst=None):
    """Load CNA data from a file and return a list of Mutation tuples.
 
    Arguments:
    cna_file -- path to TSV file containing CNAs where the first column of each line is a sample ID
                and subsequent columns contain gene names followed by "(A)" or "(D)" indicating an
                ammplification or deletion in that gene for the sample. Lines starting with '#'
                will be ignored.
    gene_wlist -- whitelist of allowed genes (default None). Genes not in this list will be ignored.
                  If None, all mutated genes will be included.
    sample_wlist -- whitelist of allowed samples (default None). Samples not in this list will be
                    ignored.  If None, all samples will be included.

    """
    arrs = [l.rstrip().split("\t") for l in open(cna_file) if not l.startswith("#")]
    return [Mutation(arr[0], cna.split("(")[0], get_mut_type(cna))
            for arr in arrs if include(arr[0], sample_wlst)
            for cna in arr[1:] if include(cna.split("(")[0], gene_wlst)]

def get_mut_type(cna):
    if cna.endswith("(A)"): return AMP
    elif cna.endswith("(D)"): return DEL
    else: raise ValueError("Unknown CNA type in '%s'", cna)

def load_oncodrive_data(fm_scores, cis_amp_scores, cis_del_scores):
    print "* Loading oncodrive data..."
    # Create defaultdicts to hold the fm and cis scores
    one = lambda: 1
    gene2fm = defaultdict(one)
    gene2cis_amp, gene2cis_del = defaultdict(one), defaultdict(one)
    
    # Load fm scores (pvals, not z-scores)
    arrs    = [ l.rstrip().split("\t") for l in open(fm_scores)
                if not l.startswith("#") ]
    gene2fm.update([(arr[1], float(arr[2])) for arr in arrs
                    if arr[2] != "" and arr[2] != "-0" and arr[2] != "-"])
    print "\tFM genes:", len(gene2fm.keys())

    # Load amplifications
    arrs = [ l.rstrip().split("\t") for l in open(cis_amp_scores)
             if not l.startswith("#")]
    gene2cis_amp.update([(arr[0], float(arr[-1])) for arr in arrs])
    print "\tCIS AMP genes:", len(gene2cis_amp.keys())

    # Load deletions
    arrs = [ l.rstrip().split("\t") for l in open(cis_del_scores)
             if not l.startswith("#")]
    gene2cis_del.update([(arr[0], float(arr[-1])) for arr in arrs])
    print "\tCIS DEL genes:", len(gene2cis_del.keys())
    
    # Merge data
    genes = set(gene2cis_del.keys()) | set(gene2cis_amp.keys()) | set(gene2fm.keys())
    print "\t- No. genes:", len(genes)
    gene2heat = dict()
    for g in genes:
        gene2heat[g] = {"del": gene2cis_del[g], "amp": gene2cis_amp[g],
                        "fm": gene2fm[g] }

    return gene2heat

def load_mutsig_scores( scores_file ):
    arrs = [l.rstrip().split("\t") for l in open(scores_file)
            if not l.startswith("#")]
    print "* Loading MutSig scores in", len(arrs), "genes..."
    return dict([(arr[0], {"pval": float(arr[-2]), "qval": float(arr[-1])})
                 for arr in arrs])


FDR_CT, FDR_LRT, FDR_FCPT = 12, 11, 10
music_score2name = {FDR_CT: "FDR_CT", FDR_LRT: "FDR_LRT", FDR_FCPT: "FDR_FCPT"}
def load_music_scores(scores_file):
    print "* Loading MuSiC scores using the median of the 3 q-values..."
    # Load file and tab-split lines 
    arrs = [l.rstrip().split("\t") for l in open(scores_file)
            if not l.startswith("#")]

    # Indices for the columns we may be interested in
    gene2music = dict([(arr[0], {"FDR_CT": float(arr[FDR_CT]),
                                 "FDR_FCPT": float(arr[FDR_FCPT]),
                                 "FDR_LRT":float(arr[FDR_LRT])})
                       for arr in arrs])

    # Output parsing info
    print "\t- Loaded %s genes." % len(gene2music)
    return gene2music
