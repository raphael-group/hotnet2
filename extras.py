import networkx as nx
strong_ccs = nx.strongly_connected_components

def coverage(genes, gene2mutations, n):
    samples = [ p for g in genes for p in gene2mutations[g].keys()
                if g in gene2mutations.keys() ]
    return len(samples) / n

def write_components_coverage(G, gene2mutations, samples, min_length=1):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs = [cc for cc in ccs if len(cc) >= min_length]
    coverages = [coverage(cc, gene2mutations, float(len(samples))) for cc in ccs]
    return "\n".join(["%s\\% %s" % (cov*100, " ".join(sorted(cc)))
                      for cov, cc in sorted(zip(coverages(ccs)))])
    
def top_genes(G, gene2heat, min_length=1, num=100):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs          = [cc for cc in ccs if len(cc) >= min_length]
    top          = sorted(gene2heat.items(), key=lambda (g, h): h)[-num:]
    genes_in_ccs = set([ g for cc in ccs for g in cc ])
    found        = [ g for g, h in top if g in genes_in_ccs ]
    not_found    = [ g for g, h in top if g not in genes_in_ccs ]
    print"Found %s of %s genes with highest heat" % (len(found), num)
    return found, not_found

def score_fn(mat, gene_index):
    gene2index = dict([(gene, index) for index, gene in gene_index.items()])
    def score(g1, g2):
        e1 = "%s -> %s = %s" % ( g1, g2, mat[gene2index[g1]][gene2index[g2]] )
        e2 = "%s -> %s = %s" % ( g2, g1, mat[gene2index[g2]][gene2index[g1]] )
        return e1 + "\n" + e2
    return score