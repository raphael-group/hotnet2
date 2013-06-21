# -*- coding: iso-8859-1 -*-
import networkx as nx, scipy as sp, numpy as np
strong_ccs = nx.strongly_connected_components

################################################################################
# Influence and similarity matrix functions

def induce_infmat(infmat, index2gene, genelist):		
    print "* Inducing infmat..."
    start_index = min(index2gene.keys())
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
    return M, index2gene

def heat_vec(gene2heat, index2gene):
    v = [gene2heat[gene] for _, gene in sorted(index2gene.iteritems()) if gene in gene2heat]
    return np.array(v)

def similarity_matrix(M, heat, directed=True):
    if directed:
        sim = M * heat
    else:
        M = np.minimum(M, M.transpose())            #ensure that the influence matrix is symmetric
        sim = np.empty_like(M)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                sim[i][j] = max(heat[i], heat[j]) * M[i][j]
    
    return sim

################################################################################
# Weighted graph functions

def weighted_graph(sim_mat, index2gene, delta, directed=True):
    e = zip( *sp.where( sim_mat >= delta) )
    edges = [(int(j), int(i), dict(weight=sim_mat[i,j])) for i, j in e]
    G = nx.DiGraph() if directed else nx.Graph()
    G.add_edges_from( [(index2gene[i], index2gene[j], d) for i, j, d in edges] )
    return G

def connected_components(G, min_size=1):
    ccs = strong_ccs(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs = [cc for cc in ccs if len(cc) >= min_size]
    return ccs

def component_sizes(ccs):
    size_dict = dict()
    for cc in ccs:
        if len(cc) not in size_dict:
            size_dict[len(cc)] = 0
        size_dict[len(cc)] += 1
    return size_dict
