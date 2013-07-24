import json

################################################################################
# Data loading functions

def load_index(index_file):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_ppi_edges(edge_list_file):
    """Load PPI edges from file and return as a set of 2-tuples of gene indices.
    
    Keyword arguments:
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
    
    Keyword arguments:
    gene_file -- path to file containing gene names, one per line.
    
    """
    return set([l.strip() for l in open(gene_file)])

def load_samples(sample_file):
    """Load sample IDs from a file and return as a set.
    
    Keyword arguments:
    sample_file -- path to TSV file containing sample IDs as the first column. Any other columns
                   will be ignored.
    
    """
    return set([l.rstrip().split()[0] for l in open(sample_file)])

def include_sample(sample, sample_wlst):
    return sample in sample_wlst if sample_wlst else True

def gene_filter(gene_set, gene_wlst):
    return gene_set & gene_wlst if gene_wlst else gene_set

def load_snv_data(snv_file, gene_wlst=None, sample_wlst=None):
    """Load SNV data from a file and return a dict mapping sample IDs to sets of genes mutated in the sample.
    
    Keyword arguments:
    snv_file -- path to TSV file containing SNVs where the first column of each line is a sample ID
                and subsequent columns are names of genes mutated in that sample. Lines starting
                with '#' will be ignored
    gene_wlist -- whitelist of allowed genes (default None). Genes not in this list will be ignored.
                  If None, all mutated genes will be included.
    sample_wlist -- whitelist of allowed samples (default None). Samples not in this list will be
                    ignored.  If None, all samples will be included.
    
    """
    arrs = [l.rstrip().split("\t") for l in open(snv_file) if not l.startswith("#")]
    return dict([(arr[0], gene_filter(set(arr[1:]), gene_wlst)) for arr in arrs
                 if include_sample(arr[0], sample_wlst)])

def load_cna_data(cna_file, gene_wlst=None, sample_wlst=None):
    """Load CNA data from a file and return a dict mapping sample IDs to dicts of (gene -> "amp"/"del") mutated in the sample.
    
    Keyword arguments:
    cna_file -- path to TSV file containing CNAs where the first column of each line is a sample ID
                and subsequent columns contain gene names followed by "(A)" or "(D)" indicating an
                ammplification or deletion in that gene for the sample. Lines starting with '#'
                will be ignored.
    gene_wlist -- whitelist of allowed genes (default None). Genes not in this list will be ignored.
                  If None, all mutated genes will be included.
    sample_wlist -- whitelist of allowed samples (default None). Samples not in this list will be
                    ignored.  If None, all samples will be included.
    
    """

    def get_mut_type(cna):
        if cna.endswith("(A)"): return "amp"
        elif cna.endswith("(D)"): return "del"
        else: raise ValueError("Unknown CNA type in sample %s: %s", arr[0], cna)
    
    arrs = [l.rstrip().split("\t") for l in open(cna_file) if not l.startswith("#")]
    arrs = [arr for arr in arrs if include_sample(arr[0], sample_wlst)]
    
    samples2cnas = {}
    for arr in arrs:
        genes = gene_filter(set([gene.split("(")[0] for gene in arr[1:]]), gene_wlst)
        gene2type = dict([(cna.split("(")[0], get_mut_type(cna)) for cna in arr[1:] if cna.split("(")[0] in genes])
        
        if len(gene2type) > 0:
            samples2cnas[arr[0]] = gene2type
            
    return samples2cnas

def load_oncodrive_data(fm_scores, cis_amp_scores, cis_del_scores):
    print "* Loading oncodrive data..."
    # Create defaultdicts to hold the fm and cis scores
    from collections import defaultdict
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
