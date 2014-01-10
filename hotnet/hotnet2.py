# -*- coding: iso-8859-1 -*-
import networkx as nx, scipy as sp, numpy as np
strong_ccs = nx.strongly_connected_components
from collections import defaultdict

################################################################################
# Influence and similarity matrix functions

def induce_infmat(infmat, index2gene, genelist):		
    """Create and return induced influence matrix containing influence scores only for those genes
    in the given gene list and a index to gene mapping for the returned matrix.
    
    Arguments:
    infmat -- 2D ndarray representing the full influence matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    genelist -- list of genes for which influence scores should be preserved in the returned matrix.
                This should be the genes that have heat scores.
     
    """
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
    """Create and return a linear ndarray of heat scores ordered by gene index.
    
    Arguments:
    gene2heat -- dict mapping a gene name to the heat score for that gene
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix. Heat scores will only be included in the returned
                  ndarray for entries in this dict.
    
    """
    v = [gene2heat[gene] for _, gene in sorted(index2gene.iteritems())]
    return np.array(v)

def similarity_matrix(M, heat):
    """Create and return a similarity matrix from the given influence matrix and heat vector.
    
    Arguments:
    M -- 2D ndarray representing the induced influence matrix obtained from induce_infmat
    heat -- 1D ndarray representing the heat score vector obtained from heat_vec
    
    """
    return M * heat

################################################################################
# Weighted graph functions

def weighted_graph(sim_mat, index2gene, delta):
    """Construct and return weighted graph in which nodes are labeled with gene names and edges
    between nodes have weight equal to the similarity score of the two genes.
    
    Arguments:
    sim_mat -- similarity matrix obtained from similarity_matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the similarity matrix
    delta -- weight threshold for edge inclusion. Only gene pairs with similarity >= delta will
             have edges connecting their corresponding nodes in the returned graph.
    
    """
    e = zip( *sp.where(sim_mat >= delta))
    edges = [(int(j), int(i), dict(weight=sim_mat[i,j])) for i, j in e]
    G = nx.DiGraph()
    G.add_edges_from([(index2gene[i], index2gene[j], d) for i, j, d in edges])
    return G

def connected_components(G, min_size=1):
    """Find connected components in the given graph and return as a list of lists of gene names.
    
    If the graph contains no connected components of size at least min_size, an empty list is returned.
    
    Arguments:
    G -- weighted graph in which connected components should be found
    min_size -- minimum size for connected components included in the returned component list 
    
    """
    ccs = [cc for cc in strong_ccs(G) if len(cc) >= min_size]
    return ccs

def component_sizes(ccs):
    """Return dict mapping a CC size to the number of connected components of that size.
    
    Only sizes for which there is at least one connected component of that size will have an entry
    in the returned dict. If the given component list is empty, an empty dict is returned.
    
    Arguments:
    ccs -- list of lists of gene names representing connected components
    
    """
    size_dict = defaultdict(int)
    for cc in ccs:
        size_dict[len(cc)] += 1
    return size_dict
