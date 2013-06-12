# -*- coding: iso-8859-1 -*-
import networkx as nx, scipy as sp, numpy as np
strong_ccs = nx.strongly_connected_components

################################################################################
# Data loading functions

def load_index( index_file ):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_heat( heat_file ):
    arrs  = [l.split() for l in open(heat_file)]
    return dict([(arr[0], float(arr[1])) for arr in arrs])

################################################################################
# Influence and similarity matrix functions

def score_fn( mat, gene_index ):
    gene2index = dict([(gene, index) for index, gene in gene_index.items()])
    def score(g1, g2):
        e1 = "%s -> %s = %s" % ( g1, g2, mat[gene2index[g1]][gene2index[g2]] )
        e2 = "%s -> %s = %s" % ( g2, g1, mat[gene2index[g2]][gene2index[g1]] )
        return e1 + "\n" + e2
    return score

def induce_infmat(infmat, index2gene, genelist, start_index=1):					
    print "* Inducing infmat..."												#shouldn't start_index be calculated as min(index2gene.keys())?
    # Reformat gene index
    gene2index = dict([(gene, index) for index, gene in index2gene.items()])

    # Identify genes in the given list that are also in the network
    genelist = [g for g in genelist if g in gene2index.keys()]
    indices = [gene2index[g]-start_index for g in genelist]
    print "\t- Genes in list and network:", len( indices )

    # Create an induced influence matrix
    M = np.zeros( (len(genelist), len(genelist)) )
    for i in range(len(indices)):
        M[i,] = infmat[gene2index[genelist[i]] - start_index, indices]

    # Create new gene index and score function
    index2gene = dict([(i, genelist[i]) for i in range(len(genelist))])
    score = score_fn( M, index2gene )
    return M, index2gene, score

def heat_vec(gene2heat, gene_index):
    v = [gene2heat[gene] for index, gene in sorted(gene_index.items()) ]
    return np.array(v)

def similarity_matrix(M, heat, gene_index, directed=True):
    if directed:
        sim = M * heat
    else:
        M = np.minimum(M, M.transpose())            #ensure that the influence matrix is symmetric
        sim = np.empty_like(M)                      #np.maximum
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                sim[i][j] = max(heat[i], heat[j]) * M[i][j]
    
    scores = score_fn( sim, gene_index )
    return sim, scores

################################################################################
# Weighted graph functions

def weighted_graph(sim_mat, gene_index, delta, directed=True):
    e = zip( *sp.where( sim_mat >= delta) )
    edges = [(int(j), int(i), dict(weight=sim_mat[i,j])) for i, j in e]
    G = nx.DiGraph() if directed else nx.Graph()
    G.add_edges_from( [(gene_index[i], gene_index[j], d) for i, j, d in edges] )
    return G

def write_components(G, min_length=1):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    return "\n".join([ " ".join(sorted(cc)) for cc in ccs
                      if len(cc) >= min_length])

def coverage(genes, gene2mutations, n):
    samples = [ p for g in genes for p in gene2mutations[g].keys()
                if g in gene2mutations.keys() ]
    return len(samples) / n

def write_components_coverage(G, gene2mutations, samples, min_length=1):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs = [cc for cc in ccs if len(cc) >= min_length]						#note: changed to >= (wherever applicable)
    coverages = [coverage(cc, gene2mutations, float(len(samples))) for cc in ccs]
    pairs = sorted(zip(coverages, ccs))
    return "\n".join(["%s\\% %s" % (cov*100, " ".join(sorted(cc)))
                      for cov, cc in sorted(zip(coverages(ccs)))])

def component_sizes(G, min_length=1):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    return ", ".join(map(str, sorted([len(cc) for cc in ccs if len(cc) >= min_length],	#rename min_length to min_size
                                     reverse=True)))

def top_genes(G, gene2heat, min_length=1, num=100):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs          = [cc for cc in ccs if len(cc) >= min_length]
    top          = sorted(gene2heat.items(), key=lambda (g, h): h)[-num:]
    genes_in_ccs = set([ g for cc in ccs for g in cc ])
    found        = [ g for g, h in top if g in genes_in_ccs ]
    not_found    = [ g for g, h in top if g not in genes_in_ccs ]
    print"Found %s of %s genes with highest heat" % (len(found), num)
    return found, not_found
