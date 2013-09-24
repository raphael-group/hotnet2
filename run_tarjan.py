import scipy.io
import hnio
import hotnet2 as hn
from tarjan import *

def run():
    infmat = scipy.io.loadmat('/data/compbio/datasets/HeatKernels/pagerank/IREFINDEX/9.0/iref_ppr_0.55.mat')['PPR']  
    infmat_index = hnio.load_index('/data/compbio/datasets/HeatKernels/pagerank/IREFINDEX/9.0/iref_index_genes')
    heat, _ = hnio.load_heat_json('/research/compbio/users/jeldridg/PanCancer/HeatFiles/pancan-freq-filtered.json')
  
    # compute similarity matrix and extract connected components
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h, True)
    G = hn.weighted_graph(sim, gene_index, 0, True)
    edges = [e(edge[0], edge[1], edge[2]['weight']) for edge in G.edges(data=True)]
    print "Finished building edges"
    edges.sort(key=sort_key)
    print "Finished sorting edges"
    T = tarjan(G, edges, 0)
    write(T, "iRef-FREQ-OUT")
    draw(T)


if __name__ == "__main__": 
    run()
