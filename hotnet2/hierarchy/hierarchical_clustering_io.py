########################################################################
#
#   Data input/output functions
#
########################################################################

def choose_labels(infmat_choice,heatvec_choice,data_choice=-1,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    if infmat_choice==0:
        infmat_label = 'HINT+HI2012'
    elif infmat_choice==1:
        infmat_label = 'iRefIndex'
    elif infmat_choice==2:
        infmat_label = 'Multinet'

    if heatvec_choice==0:
        heatvec_label = 'MutationFrequency'
    elif heatvec_choice==1:
        heatvec_label = 'MutSigCV'

    if data_choice==0:
        data_label = 'real'
    elif data_choice==1:
        data_label = 'uniform'
    elif data_choice==2:
        if permutation_choice>-1:
            data_label = 'permuted_%d' % (permutation_choice+1)
        else:
            raise Exception('The permutation number for the permuted heat vector is unspecified or invalid.')
    elif data_choice==3:
        if permutation_choice>-1:
            data_label = 'implanted_%d_%d_%d_%d' % (permutation_choice+1,implant_selection_choice+1,implant_number_choice,implant_size_choice)
        else:
            raise Exception('The permutation number for the permuted heat vector is unspecified or invalid.')

    if data_choice==-1:
        return infmat_label, heatvec_label
    else:
        return infmat_label, heatvec_label, data_label

def import_data(infmat_choice,heatvec_choice):

    import numpy as np
    import scipy
    import scipy.io

    if infmat_choice==0:
        infmat_filename = '/data/compbio/datasets/HeatKernels/pagerank/IntHint/binary+complex+hi2012/inthint_ppr_0.60.mat'
        infmat_type = 'PPR'
        geneindex_filename = '/data/compbio/datasets/HeatKernels/pagerank/IntHint/binary+complex+hi2012/inthint_index_genes'
    elif infmat_choice==1:
        infmat_filename = '/data/compbio/datasets/HeatKernels/pagerank/IREFINDEX/9.0/iref_ppr_0.55.mat'
        infmat_type = 'PPR'
        geneindex_filename = '/data/compbio/datasets/HeatKernels/pagerank/IREFINDEX/9.0/iref_index_genes'
    elif infmat_choice==2:
        infmat_filename = '/data/compbio/datasets/HeatKernels/pagerank/MultiNet/1.0/multinet_ppr_0.50.mat'
        infmat_type = 'PPR'
        geneindex_filename = '/data/compbio/datasets/HeatKernels/pagerank/MultiNet/1.0/multinet_index_genes'

    if heatvec_choice==0:
        heatvec_filename = '/research/compbio/projects/HotNet2/ng/data/scores/pan12.gene2freq.txt'
    elif heatvec_choice==1:
        heatvec_filename = '/research/compbio/projects/HotNet2/ng/data/scores/pan12.gene2mutsig.txt'

    indices_to_genes = {}
    genes_to_indices = {}
    with open(geneindex_filename) as geneindex_file:
        for line in geneindex_file:
            (index, gene) = line.split()[0:2]
            indices_to_genes[int(index)-1] = gene
            genes_to_indices[gene] = int(index)-1

    genes_to_heats = {}
    with open(heatvec_filename) as heatvec_file:
        for line in heatvec_file:
            (gene, heat) = line.split()
            genes_to_heats[gene] = float(heat)

    common_genes = sorted(list(set.intersection(set(genes_to_indices.keys()), set(genes_to_heats.keys()))))

    indices = np.array([genes_to_indices[gene] for gene in common_genes], dtype = np.int)
    infmat = np.array(scipy.io.loadmat(infmat_filename)[infmat_type], dtype=np.float)
    infmat = infmat[np.ix_(indices,indices)]
    heats = np.array([genes_to_heats[gene] for gene in common_genes], dtype = np.float)

    return infmat, heats, common_genes

def import_tree(infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    path = '/data/compbio/users/mreyna/decomposition-trees/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    filename = path + 'tree'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'

    with open(filename,'rb') as pickle_file:
        T = pickle.load(pickle_file)

    return T

def export_tree(T,infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    path = '/data/compbio/users/mreyna/decomposition-trees/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    filename = path + 'tree'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'

    with open(filename,'wb') as pickle_file:
        pickle.dump(T,pickle_file,protocol=-1)

def import_weight_cluster(infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    path = '/data/compbio/users/mreyna/weights-clusters/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    weights_filename = path + 'weight'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'
    clusters_filename = path + 'cluster'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'

    with open(weights_filename,'rb') as pickle_file:
        weights = pickle.load(pickle_file)
    with open(clusters_filename,'rb') as pickle_file:
        clusters = pickle.load(pickle_file)

    return weights, clusters

def export_weight_cluster(weights,clusters,infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    path = '/data/compbio/users/mreyna/weights-clusters/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    weights_filename = path + 'weight'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'
    clusters_filename = path + 'cluster'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.pickle'

    with open(weights_filename,'wb') as pickle_file:
        pickle.dump(weights,pickle_file,protocol=-1)
    with open(clusters_filename,'wb') as pickle_file:
        pickle.dump(clusters,pickle_file,protocol=-1)

def import_linkage(infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    import json

    path = '/data/compbio/users/mreyna/linkage-matrices/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    filename = path + 'linkage'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.json'

    with open(filename,'r') as json_file:
        contents = json.load(json_file)

    Z = contents["Z"]
    labels = contents["labels"]

    return Z, labels

def export_linkage(Z,labels,infmat_choice,heatvec_choice,data_choice,permutation_choice=-1,implant_selection_choice=-1,implant_number_choice=-1,implant_size_choice=-1):

    import json

    path = '/data/compbio/users/mreyna/linkage-matrices/'
    infmat_label, heatvec_label, data_label = choose_labels(infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
    filename = path + 'linkage'+'_'+infmat_label+'_'+heatvec_label+'_'+data_label+'.json'

    with open(filename,'w') as json_file:
        json.dump({"Z":Z,"labels":labels},json_file)

########################################################################
#
#   Output conversion functions
#
########################################################################

# The next function converts a tree from our specialized format to a linkage
# matrix format with linkage matrix Z and leaf node list V.

def linkage(T):

    condensations = sorted([v for v in T if len(v)==2])
    increasing = condensations[0][0]<T[condensations[0]][0]
    inner_nodes = sorted(list(set(T.values())), key=lambda e: e[0], reverse=not increasing)
    base = condensations[0][0]

    V = [v[1] for v in condensations]
    L = {v:i for i,v in enumerate(condensations)}
    Z = []
    k = len(condensations)

    for w in inner_nodes:

        children = [v for v in condensations if T[v]==w]

        x = children[0]
        for y in children[1:]:
            z = tuple([w[0]]+sorted(x[1:]+y[1:]))
            Z.append([L[x],L[y],w[0],len(z[1:])])
            L[z] = k
            x = z
            k += 1

        for v in children:
            condensations.remove(v)
        condensations.append(w)

    # Using the label ordering provided by SciPy's dendrogram function,
    # reorder the leaf nodes in Z and V.
    import  scipy.cluster, sys
    sys.setrecursionlimit(100000)

    Y = [[a,b,base-c,d] for (a,b,c,d) in Z]
    R = scipy.cluster.hierarchy.dendrogram(Y,labels=V,orientation='right', no_plot=True)
    W = list(reversed(R["ivl"]))

    n = len(W)
    for row in Y:
        if row[0] < n:
            row[0] = W.index(V[row[0]])
        if row[1] < n:
            row[1] = W.index(V[row[1]])

    return Y,W

# The next function converts a tree from our specialized format to the
# standard Newick format.

def newick(T):

    condensations = [v for v in T if len(v)==2]
    increasing = condensations[0][0]<T[condensations[0]][0]
    inner_nodes = sorted(list(set(T.values())), key=lambda e: e[0], reverse=not increasing)
    root = inner_nodes.pop()

    mapping = {v:str(v[1])+':'+str(abs(T[v][0]-v[0])) for v in condensations}

    for w in inner_nodes:
        children = [v for v in condensations if T[v]==w]
        mapping[w] = '('+','.join([mapping[v] for v in children])+'):'+str(abs(T[w][0]-w[0]))
        for v in children:
            condensations.remove(v)
        condensations.append(w)

    children = [v for v in condensations if T[v]==root]
    return '('+','.join([mapping[v] for v in children])+');'

########################################################################
#
#   Progress/debugging input/output functions
#
########################################################################

def progress(message):

    import sys

    try:
        previous_length = progress.length
    except:
        previous_length = 0

    sys.stdout.write("\r"+" "*(previous_length))
    sys.stdout.flush()
    sys.stdout.write("\r"+str(message))
    sys.stdout.flush()

    progress.length = len(str(message))

def prettify(X):

    # This function formats various data types in an easier-to-read way.

    if type(X) in [dict]:
        string = '\n'.join([''.join(['    ',str(x),' : ',str(f(x))]) for x in X])
    elif type(X) in [list, tuple]:
        string = '\n'.join([''.join(['    ',str(i),' : ',str(x)]) for i,x in enumerate(X)])
    else:
        string = str(X)

    return string

if __name__ == "__main__":

    from hierarchical_clustering import *
    import sys

    # Options for influence matrix:
    # HINT+HI2012: 0
    # iRefIndex: 1
    # Multinet: 2

    # Options for mutation frequency:
    # Mutation frequency: 0
    # MutSigCV: 1

    # Options for variations on heat vector data:
    # Real: 0
    # Uniform: 1
    # Randomly permuted: 2
    # Simulated: 3

    # Options for permutations:
    # number of permutations of heat vector: any nonnegative integer; sets random seed

    ####################################################################

    # Import command line arguments.

    if len(sys.argv)==4:
        infmat_choice, heatvec_choice, data_choice = [int(i) for i in sys.argv[1:]]
    elif len(sys.argv)==5:
        infmat_choice, heatvec_choice, data_choice, permutation_choice = [int(i) for i in sys.argv[1:]]
    elif len(sys.argv)==8:
        infmat_choice, heatvec_choice, data_choice, permutation_choice, implant_selection_choice, implant_number_choice, implant_size_choice = [int(i) for i in sys.argv[1:]]
    else:
        raise Exception('%d arguments specified; 3, 4, or 7 expected' % (len(sys.argv)-1))

    # Check validity of arguments.

    if infmat_choice not in [0,1,2]:
        raise Exception('%d given in first argument; 0, 1, or 2 expected for choice of influence network' % infmat_choice)
    if heatvec_choice not in [0,1]:
        raise Exception('%d given in second argument; 0 or 1 expected for choice of heat scores' % heatvec_choice)
    if data_choice not in [0,1,2,3]:
        raise Exception('%d given in third argument; 0, 1, 2, or 3 expected for choice of heat score modification' % data_choice)
    if data_choice==2:
        if len(sys.argv)!=5:
            raise Exception('Correct number of arguments not given; nonnegative integer expected for choice of heat score permutation')
        else:
            if permutation_choice<0:
                raise Exception('%d given in fourth argument; nonnegative integer expected for choice of heat score permutation' % permutation_choice)
    if data_choice==3:
        if len(sys.argv)!=8:
            raise Exception('Correct number of arguments not given; nonnegative integer expected for choice of heat score permutation')
        else:
            if permutation_choice<0:
                raise Exception('%d given in fourth argument; nonnegative integer expected for choice of heat score permutation' % permutation_choice)

    ####################################################################

    # Import real data.

    infmat, heats, genes = import_data(infmat_choice,heatvec_choice)

    # (Possibly) modify heat scores.

    if data_choice==0:

        sim = infmat*heats
        T = HD(genes,sim)
        Z, labels = linkage(T)

    elif data_choice==1:

        uniform_heats = np.mean(heats)*np.ones(len(heats))
        sim = infmat*uniform_heats
        T = HD(genes,sim)
        Z, labels = linkage(T)

    elif data_choice==2:

        np.random.seed(seed=permutation_choice)
        random_heats = np.random.permutation(heats)
        sim = infmat*random_heats
        T = HD(genes,sim)
        Z, labels = linkage(T)

    elif data_choice==3:

        S = import_tree(infmat_choice,heatvec_choice,1)
        clusters, delta = cluster_cutoffs(S,implant_selection_choice,implant_number_choice,implant_size_choice)
        S = None

        hot_indices = [genes.index(v) for w in clusters for v in w]
        cold_indices = [i for i in xrange(len(heats)) if i not in hot_indices]

        if permutation_choice==0:
            infmat_label, heatvec_label = choose_labels(infmat_choice,heatvec_choice)
            output = []
            output.append('Interaction network:')
            output.append('    '+infmat_label)
            output.append('Heat scores:')
            output.append('    '+heatvec_label)
            output.append('Selection critera:')
            if implant_selection_choice==0:
                output.append('    '+'%d cluster(s) of size %d or larger' %(implant_number_choice,implant_size_choice))
            elif implant_selection_choice==1:
                output.append('    '+'%d or more vertices in cluster(s) of size %d or larger' %(implant_number_choice,implant_size_choice))
            output.append('Largest delta satisfying criteria:')
            output.append('    %.4e' % delta)
            output.append('Implanted cluster(s):')
            output.append(prettify(clusters))
            output.append('Implanted cluster size(s):')
            output.append('    '+str(map(len,clusters)))
            output.append('Sum of implanted cluster size(s):')
            output.append('    '+str(sum(map(len,clusters))))
            print '\n'.join(output)
   
        sorted_heats = np.sort(heats)
        implanted_heats = np.empty(len(heats))
        np.random.seed(seed=permutation_choice)
        implanted_heats[hot_indices] = np.random.permutation(sorted_heats[-len(hot_indices):])
        implanted_heats[cold_indices] = np.random.permutation(sorted_heats[:-len(hot_indices)])        

        sim = infmat*implanted_heats
        T = HD(genes,sim)
        Z, labels = linkage(T)      

    # Compute hierarchical decomposition tree and linkage matrix and labels.

    # Store results.

    if data_choice==0 or data_choice==1:
        export_tree(T,infmat_choice,heatvec_choice,data_choice)
        export_linkage(Z,labels,infmat_choice,heatvec_choice,data_choice)
    elif data_choice==2:
        export_tree(T,infmat_choice,heatvec_choice,data_choice,permutation_choice)
        export_linkage(Z,labels,infmat_choice,heatvec_choice,data_choice,permutation_choice)
    elif data_choice==3:
        export_tree(T,infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
        export_linkage(Z,labels,infmat_choice,heatvec_choice,data_choice,permutation_choice,implant_selection_choice,implant_number_choice,implant_size_choice)
