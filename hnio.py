import json

################################################################################
# Data loading functions

def load_index(index_file):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_heat_json(heat_file):
    with open(heat_file) as f:
        blob = json.load(f)
        return blob["heat"], blob["parameters"]

def load_heat_tsv(heat_file):
    arrs  = [l.split() for l in open(heat_file)]
    return dict([(arr[0], float(arr[1])) for arr in arrs])

def load_gene_list(gene_list_file):
    return set([l.strip() for l in open(gene_list_file)])

def load_mutation_data(snv_file, cna_file, sample_file, gene_file):
    print "* Loading mutation data..."
    # Load whitelist
    gene_wlst = set([ l.rstrip() for l in open(gene_file) ])
    sample_wlst = set([ l.rstrip().split()[0] for l in open(sample_file) ])
    samples_w_tys = [ tuple( l.rstrip().split()) for l in open(sample_file) ]

    # Load SNVs
    print "\tLoading SNVs..."
    arrs = [l.rstrip().split("\t") for l in open(snv_file)
            if not l.startswith("#")]
    sample2snvs = dict([(arr[0], set(arr[1:]) & gene_wlst) for arr in arrs
                        if arr[0] in sample_wlst])
    samples = set( sample2snvs.keys() )
    genes = set([ g for _, snvs in sample2snvs.items() for g in snvs])
    genes = genes & gene_wlst
    gene2heat = dict([(g, dict()) for g in genes])
    for sample, snvs in sample2snvs.items():
        for g in snvs:
            gene2heat[g][sample] = [ "snv" ]

    # Load CNAs
    print "\tLoading CNAs..."
    arrs = [l.rstrip().split("\t") for l in open(cna_file)
            if not l.startswith("#")]
    arrs = [arr for arr in arrs if arr[0] in sample_wlst]
    
    new_genes = set([ g.split("(")[0] for a in arrs for g in a[1:]])
    new_genes = new_genes & gene_wlst
    for g in set( new_genes ) - set( gene2heat.keys() ): gene2heat[g] = dict()
    genes = set(gene2heat.keys())
    for arr in arrs:
        sample, mutations = arr[0], set([g for g in arr[1:]
                                          if g.split("(")[0] in gene2heat.keys() ])
        samples.add( sample )
        for g in mutations:
            
            if g.endswith("(A)"): mut_type = "amp"
            elif g.endswith("(D)"): mut_type = "del"
            gene_name = g.split("(")[0]
                
            if sample in gene2heat[gene_name].keys():
                gene2heat[gene_name][sample].append( mut_type )
            else:
                gene2heat[gene_name][sample] = [ mut_type ]
    genes = set(gene2heat.keys())

    print "\tLoaded mutation data for", len(genes),
    print "genes in %s samples" % len(samples)

    sample_w_tys = [(s, ty) for s, ty in samples_w_tys if s in samples]     #MAX: what's the deal with this line?

    return (genes, samples, gene2heat), samples_w_tys

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
