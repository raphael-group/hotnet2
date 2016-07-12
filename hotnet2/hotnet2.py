# -*- coding: iso-8859-1 -*-
from collections import defaultdict
import networkx as nx, numpy as np, scipy as sp

try:
    import c_routines
    fast_similarity_matrix = True
except ImportError:
    print("WARNING: Could not import either C module; "
          "falling back to NumPy for similarity matrix creation.")
    fast_similarity_matrix = False

################################################################################
# Influence and similarity matrix functions

def similarity_matrix(infmat, index2gene, gene2heat, directed=True, verbose=0):
    """Create and return a similarity matrix and index to gene mapping for the given influence
    matrix and heat. Only genes with heat that are in the network will be included in the returned
    similarity matrix and index to gene mapping.

    Arguments:
    infmat -- 2D ndarray representing the full influence matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    gene2heat -- dict mapping a gene name to the heat score for that gene
    directed -- if True, sim[i][j] = inf(i,j)*heat[i] and sim[i][j] != sim[j][i]
                if False, sim[i][j] = min(inf(i,j), inf(j,i))*max(heat(i), heat(j))

    """
    start_index = min(index2gene.keys())
    gene2index = dict((gene, index) for index, gene in index2gene.iteritems())

    # Identify genes in the given list that are also in the network
    genelist = sorted(set(gene2heat.keys()).intersection(gene2index.keys()))
    index2gene = dict(enumerate(genelist))
    if verbose > 4:
        print "\t- Genes in similarity matrix:", len(genelist)

    infmat = np.asarray(infmat, dtype=np.float64)
    h = np.array([gene2heat[g] for g in genelist], dtype=np.float64)
    indices = np.array([gene2index[g]-start_index for g in genelist], dtype=np.int)
    m = np.shape(infmat)[0]
    n = np.shape(h)[0]

    if fast_similarity_matrix:
        if directed:
            sim = c_routines.compute_sim(infmat, h, indices, m, n)
        else:
            sim = c_routines.compute_sim_classic(infmat, h, indices, m, n)
    else:
        M = infmat[np.ix_(indices, indices)]
        if directed:
            sim = M * h
        else:
            M = np.minimum(M, M.transpose())  # Ensure that the influence matrix is symmetric
            sim = np.empty_like(M)
            for i in range(n):
                for j in range(i, n):
                    sim[i, j] = max(h[i], h[j]) * M[i, j]
                    sim[j, i] = sim[i, j]

    return sim, index2gene

################################################################################
# Weighted graph functions

def weighted_graph(sim_mat, index2gene, delta, directed=True):
    """Construct and return weighted graph in which nodes are labeled with gene names and edges
    between nodes have weight equal to the similarity score of the two genes.

    Arguments:
    sim_mat -- similarity matrix obtained from similarity_matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the similarity matrix
    delta -- weight threshold for edge inclusion. Only gene pairs with similarity >= delta will
             have edges connecting their corresponding nodes in the returned graph.
    directed -- whether or not the resulting graph should be directed. If true, a networkx Digraph
                instance is returned. If false, a networkx Graph instance is returned.

    """
    e = zip( *sp.where(sim_mat >= delta))
    edges = [(int(j), int(i), dict(weight=sim_mat[i,j])) for i, j in e]
    G = nx.DiGraph() if directed else nx.Graph()
    G.add_edges_from([(index2gene[i], index2gene[j], d) for i, j, d in edges])
    return G

def connected_components(G, min_size=1):
    """Find connected components in the given graph and return as a list of lists of gene names.

    If the graph contains no connected components of size at least min_size, an empty list is returned.

    Arguments:
    G -- weighted graph in which connected components should be found
    min_size -- minimum size for connected components included in the returned component list

    """
    ccs = nx.strongly_connected_components(G) if isinstance(G, nx.DiGraph) else nx.connected_components(G)
    ccs = [cc for cc in ccs if len(cc) >= min_size]
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
