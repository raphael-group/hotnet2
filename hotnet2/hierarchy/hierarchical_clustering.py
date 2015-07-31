#!/usr/bin/python

###############################################################################
#
#   Setup
#
###############################################################################

import math
import numpy as np

try:
    try:
        from .. import fortran_routines
        available_routines = 1
    except ImportError:
        import fortran_routines
        available_routines = 1
except ImportError:
    available_routines = 0

###############################################################################
#
#   Tarjan's hierarchical decomposition algorithm
#
###############################################################################

#
# HD(V,A,increasing=False)
#
#   This function finds the decomposition of a graph into a hierarchy of
#   strongly connected components following the algorithm of Tarjan (1983).
#   It is a wrapper for tarjan_HD, which actually performs the decomposition.
#
#   Inputs:
#       V: list of vertices of the graph
#       A: weighted adjacency matrix representing the edges of the graph
#       increasing: order for adding edges according to weight
#
#   Output:
#       T: hierarchical decomposition tree
#

def HD(V,A,increasing=False):

    SCCs = strongly_connected_components(A)

    if len(SCCs)>1:
        largest_component_size = max(map(len,SCCs))
        number_vertices = len(V)
        for component in SCCs:
            if len(component)==largest_component_size:
                break
        A = A[np.ix_(component,component)]
        V = slice_vertices(V,component)
        print 'The graph has %d strongly connected components, but the hierarchical decomposition requires a strongly connected graph.  The decomposition tree considers a subgraph defined by a largest strongly connected component of the original graph %d of its %d vertices.' % (len(SCCs),largest_component_size,number_vertices)

    if increasing:
        base = 0.0
        W = [tuple([base,v]) for v in V]
        T,root = tarjan_HD(W,A,{},0)
    elif not increasing:
        base = np.max(A)
        W = [tuple([base,v]) for v in V]
        S,root = tarjan_HD(W,reverse_matrix(A),{},0)
        T = reverse_tree(W,A,S)

    return T

#
# tarjan_HD(V,A,T,i)
#
#   This function implements the hierarchical decomposition algorithm
#   described in Tarjan (1983).  Whenever possible, we reuse the notation and
#   labels, i.e., Case 1, Case 2, Case 2a, and Case 2b, used in the paper.
#
#   Inputs:
#       V: list of vertices of the graph
#       A: weighted adjacency matrix representing the edges of the graph
#       T: hierarchical decomposition tree
#       i: weight index; see Tarjan (1983) for description
#
#   Outputs:
#       T: hierarchical decomposition tree
#       root: root of the tree
#

def tarjan_HD(V,A,T,i):

    A_sort = sort(A)
    m = len(A_sort)-1
    r = m-i

    if r==1:
        # Case 1
        root = form_edge(A_sort[-1],V)
        S = {v:root for v in V}
        T.update(S)
        return T,root

    else:
        # Case 2
        j = int(math.ceil(0.5*float(i+m)))
        weight_i = A_sort[i]
        weight_j = A_sort[j]
        A_j = remove_edges(A,weight_j)
        SCCs = strongly_connected_components(A_j)

        if len(SCCs)==1:
            # Case 2a
            return tarjan_HD(V,A_j,T,i)

        else:
            # Case 2b
            W = []
            for SCC in SCCs:
                if len(SCC)>1:
                    B = slice_array(A_j,SCC,SCC)
                    k = subproblem_index(B,weight_i)
                    X = slice_vertices(V,SCC)
                    S,root = tarjan_HD(X,B,{},k)
                    T.update(S)
                    W.append(root)
                else:
                    W.extend(slice_vertices(V,SCC))

            B = condense_graph(A,SCCs)
            k = subproblem_index(B,weight_j)
            return tarjan_HD(W,B,T,k)

#
# clustering(T)
#
#   This function finds the clusters of vertices given in the hierarchical
#   decomposition tree T.  Each vertex of the original graph forms its own
#   cluster at the leaf nodes of the tree while the same vertices form a single
#   cluster at the root of the tree.
#
#   Input:
#       T: hierarchical decomposition tree
#
#   Outputs:
#       weights: weights of the inner nodes of the tree
#       clusters: clusters of vertices corresponding to each weight
#

def clustering(T):

    condensations = [v for v in T if len(v)==2]
    increasing = condensations[0][0]<T[condensations[0]][0]
    clusters = [sorted([[v[1]] for v in condensations])]
    inner_nodes = sorted(list(set(T.values())), key=lambda e: e[0], reverse=not increasing)
    weights = sorted(list(set([e[0] for e in inner_nodes])), reverse=not increasing)

    i = 0
    n = len(inner_nodes)

    for delta in weights:

        while i<n and inner_nodes[i][0]==delta:
            children = [v for v in condensations if T[v]==inner_nodes[i]]
            for v in children:
                condensations.remove(v)
            condensations.append(inner_nodes[i])
            i += 1

        leaf_nodes = [list(v[1:]) for v in condensations]
        clusters.append(sorted(leaf_nodes))

    if increasing:
        weights.insert(0,0.0)
    elif not increasing:
        weights.append(0.0)

    return weights,clusters

#
# delta_clustering(T,delta)
#
#   This function finds the clusters of vertices given in the hierarchical
#   decomposition tree T for a particular delta.  Compare this to cluster,
#   which finds the clusters for all values of delta as well as the values of
#   delta that condense the components.
#
#   Input:
#       T: hierarchical decomposition tree
#
#   Output:
#       clusters: clusters of vertices corresponding to delta
#

def delta_clustering(T,delta):

    condensations = [v for v in T if len(v)==2]
    increasing = condensations[0][0]<T[condensations[0]][0]
    inner_nodes = sorted(list(set(T.values())), key=lambda e: e[0], reverse=not increasing)
    weights = sorted(list(set([e[0] for e in inner_nodes])), reverse=not increasing)

    i = 0
    n = len(inner_nodes)

    for weight in weights:

        if increasing and delta<weight:
            break
        elif not increasing and delta>weight:
            break

        while i<n and inner_nodes[i][0]==weight:
            children = [v for v in condensations if T[v]==inner_nodes[i]]
            for v in children:
                condensations.remove(v)
            condensations.append(inner_nodes[i])
            i += 1

    clusters = [list(v[1:]) for v in condensations]

    return sorted(clusters)

#
# cluster_cutoffs(T)
#
#   This function finds the clusters of vertices given in the hierarchical
#   decomposition tree T with a minimum number and minimum size.
#
#   Input:
#       T: hierarchical decomposition tree
#       selection_criterion: whether minimum number corresponds to number of
#           clusters of satisfying criteria or number of vertices satisfying
#           criteria
#       minimum_number: minimum number of clusters of given size
#       minimum_size: minimum size of clusters
#       truncated: whether or not to include singletons (default is not)
#
#   Outputs:
#       clusters: first partition of vertices satisfying criteria
#       delta: first weight at which criteria are satisfied
#

def cluster_cutoffs(T,selection_criterion,minimum_number,minimum_size,truncated=True):

    condensations = [v for v in T if len(v)==2]
    increasing = condensations[0][0]<T[condensations[0]][0]
    inner_nodes = sorted(list(set(T.values())), key=lambda e: e[0], reverse=not increasing)
    weights = sorted(list(set([e[0] for e in inner_nodes])), reverse=not increasing)

    i = 0
    n = len(inner_nodes)

    for delta in weights:

        while i<n and inner_nodes[i][0]==delta:
            children = [v for v in condensations if T[v]==inner_nodes[i]]
            for v in children:
                condensations.remove(v)
            condensations.append(inner_nodes[i])
            i += 1

        counts = [len(v[1:]) for v in condensations if len(v[1:])>=minimum_size]

        if selection_criterion==0:
            if len(counts)>=minimum_number:
                if truncated:
                    clusters = [v[1:] for v in condensations if len(v[1:])>=minimum_size]
                if not truncated:
                    clusters = [v[1:] for v in condensations]
                break
        elif selection_criterion==1:
            if sum(counts)>=minimum_number:
                if truncated:
                    clusters = [v[1:] for v in condensations if len(v[1:])>=minimum_size]
                if not truncated:
                    clusters = [v[1:] for v in condensations]
                break

    return clusters, delta

###############################################################################
#
#   Helper functions
#
###############################################################################

def closest(Y,x):

    m = len(Y)//2
    n = len(Y)//2
    r = len(Y)-1

    while m>1:
        m = int(round(m/2.0))
        if x<=Y[n]:
            n = max(n-m,0)
        else:
            n = min(n+m,r)

    p = max(n-1,0)
    q = min(n+1,r)

    if Y[p]<=x<=Y[n]:
        if abs(x-Y[p])<=abs(x-Y[n]):
            index = p
        else:
            index = n
    else:
        if abs(x-Y[n])<=abs(x-Y[q]):
            index = n
        else:
            index = q

    return index

def closest_reversed(Y,x):

    return len(Y)-1-closest(Y[::-1],x)

def condense_graph(A,SCCs):

    if available_routines==1:

        vertices = np.array([vertex for component in SCCs for vertex in component],dtype=np.int)
        sizes = [len(component) for component in SCCs]
        indices = np.array([sum(sizes[:i]) for i in xrange(len(SCCs)+1)],dtype=np.int)
        return fortran_routines.condense_graph(A,vertices+1,indices+1)

    else:

        n = len(SCCs)
        B = np.zeros((n,n),dtype=np.float)

        for i in xrange(n):
             for j in xrange(n):
                if i!=j:
                    C = slice_array(A,SCCs[j],SCCs[i])
                    D = np.nonzero(C)
                    if np.size(D)>0:
                        B[i,j] = np.min(C[D])

        return B

def form_edge(weight,V):

    W = []
    for v in V:
        W.extend(v[1:])
    return tuple([weight]+sorted(W))

def remove_edges(A,weight):

    if available_routines==1:

        return fortran_routines.remove_edges(A,weight)

    else:

        B = A.copy()
        B[B>weight] = 0
        return B

def slice_vertices(vertices,indices):

    return [vertices[index] for index in indices]

# The next function is a hack to reverse the order of the edge weights in the
# adjacency matrix: the largest weight is mapped to the smallest weight, the
# second largest to the second smallest, etc.  The reverse_tree function
# effectively reverses the reverse_matrix function.

def reverse_matrix(A):

    shift = int(math.ceil(np.max(A)))+1
    B = -A+shift
    B[B==shift] = 0
    return B

# The next function is a hack to reverse the order of the edge weights in the
# hierarchical.  Due to floating-point precision errors, it is also necessary
# to match each reversed edge weight with its original edge weight in A.  The
# reverse_tree function effectively reverses the reverse_matrix function.

def reverse_tree(V,A,S):

    weights = np.unique(A)
    interior_nodes = set(S.values())
    shift = int(math.ceil(weights[-1]))+1
    mapping = {v: v for v in V}

    for e in interior_nodes:
        u = -e[0]+shift
        w = weights[closest(weights,u)]
        mapping[e] = tuple([w]+list(e[1:]))

    return {mapping[v]: mapping[S[v]] for v in S}

def slice_array(A,rows,columns):

    if available_routines==1:

        return fortran_routines.slice_array(A,np.array(columns,dtype=np.int)+1,np.array(rows,dtype=np.int)+1)

    else:

        return A[np.ix_(rows,columns)]

def sort(A):

    return np.unique(A)

def strongly_connected_components(A):

    if available_routines==1:

        components = fortran_routines.strongly_connected_components(A)
        return [np.where(components==i+1)[0].tolist() for i in xrange(np.max(components))]

    else:

        components = strongly_connected_components_from_matrix(A)
        return sorted([sorted(c) for c in components])

# I adapted the next function from NetworkX's documentation to use adjacency
# matrices instead of adjacency lists.

def strongly_connected_components_from_matrix(A):

    n=len(A)
    preorder={}
    lowlink={}
    scc_found={}
    scc_queue = []
    i=0
    for source in xrange(n):
        if source not in scc_found:
            queue=[source]
            while queue:
                v=queue[-1]
                if v not in preorder:
                    i=i+1
                    preorder[v]=i
                done=True
                v_nbrs=[w for w in xrange(n) if A[v,w]!=0]
                for w in v_nbrs:
                    if w not in preorder:
                        queue.append(w)
                        done=False
                        break
                if done==True:
                    lowlink[v]=preorder[v]
                    for w in v_nbrs:
                        if w not in scc_found:
                            if preorder[w]>preorder[v]:
                                lowlink[v]=min([lowlink[v],lowlink[w]])
                            else:
                                lowlink[v]=min([lowlink[v],preorder[w]])
                    queue.pop()
                    if lowlink[v]==preorder[v]:
                        scc_found[v]=True
                        scc=[v]
                        while scc_queue and preorder[scc_queue[-1]]>preorder[v]:
                            k=scc_queue.pop()
                            scc_found[k]=True
                            scc.append(k)
                        yield scc
                    else:
                        scc_queue.append(v)

def subproblem_index(B,A_weight):

    C = np.unique(B)
    index = closest(C,A_weight)
    maximum_index = len(C)-1

    while C[index]<A_weight and index<maximum_index:
        index += 1
    while C[index]>A_weight and index>0:
        index -= 1
    return index

if __name__ == "__main__":

    ###########################################################################
    #
    #   Testing functions
    #
    ###########################################################################

    def edges_to_matrix(E):

        V = edges_to_vertices(E)
        n = len(V)
        indices = {V[i]:i for i in xrange(n)}

        A = np.zeros((n,n),dtype=np.float)
        for (source,target,weight) in E:
            A[indices[source],indices[target]] = weight

        return V,A

    def edges_to_vertices(E):

        V = set()

        for (source,target,weight) in E:
            V.add(source)
            V.add(target)

        return sorted(list(V))

    # The next function computes the hierarchical decomposition naively, i.e.,
    # by adding one edge at a time. The disadvantage of the naive approach
    # compared to Tarjan's algorithm is the need to compute strongly connected
    # components many more times.

    def HD_naive(E,increasing):

        import networkx as nx

        vertices = edges_to_vertices(E)
        edges = sorted(E, key = lambda e:e[2], reverse = not increasing)
        weights = sorted(list(set(e[2] for e in edges)), reverse = not increasing)

        # Add all of the vertices to the graph.

        G = nx.DiGraph()
        G.add_nodes_from(vertices)

        T = {}
        if increasing:
            roots = {v:tuple([0.0,v]) for v in vertices}
        elif not increasing:
            roots = {v:tuple([max(weights),v]) for v in vertices}

        i = 0
        m = len(edges)
        n = len(vertices)

        # Consider each weight.

        for weight in weights:

            # Add all of the edges with the same weight to the graph.

            while edges[i][2]==weight:
                G.add_edge(edges[i][0], edges[i][1], weight=edges[i][2])
                i += 1
                if i==m:
                    break

            # Find the SCCs of the resulting graph.

            components = [component for component in nx.strongly_connected_components(G)]

            # Check if any the components contracted by comparing the number of
            # components previously to the number of components currently.

            if len(components)<n:
                n = len(components)
                for component in components:

                    # Check each component for contractions by checking if
                    # every vertex in the component has the same root.  If not,
                    # a contraction occurred while adding the latest weight, so
                    # update the tree and roots.

                    first_root = roots[component[0]]
                    if not all(roots[v]==first_root for v in component):
                        root = tuple([weight]+sorted(component))
                        for v in component:
                            T[roots[v]] = root
                            roots[v] = root

        return T

    def matrix_to_edges(V,A):

        (m,n) = np.shape(A)
        o = len(V)

        if m!=n or n!=o:
            raise Exception('Number of vertices and adjacency matrix dimensions do not match.')

        return [[V[i],V[j],A[i,j]] for i in xrange(n) for j in xrange(n) if A[i,j]!=0]

    # The next function formats various data types in an arguably easier-to-
    # read way.

    def prettify(X):

        if type(X) in [dict]:
            string = ''
            for x in X:
                string += '    '+str(x)+' : '+str(X[x])+'\n'
        elif type(X) in [list,tuple]:
            string = ''
            for i,x in enumerate(X):
                string += '    '+str(i)+' : '+str(x)+'\n'
        else:
            string = '    '+str(X)+'\n'

        return string

    # The next function displays progress messages, but I have only tested it
    # on Debian and Ubuntu.

    def progress(message):

        import sys

        try:
            previous_length = len(progress.previous)
        except:
            previous_length = 0

        sys.stdout.write("\r"+" "*previous_length)
        sys.stdout.flush()
        sys.stdout.write("\r"+str(message))
        sys.stdout.flush()

        progress.previous = message

    # The next function generates weighted adjacency matrices with uniformly
    # distributed edge weights that conform to given sparsity and
    # nonuniqueness parameters.

    def random_adjacency_matrix(n,seed=np.random.randint(0,4294967295),sparsity=0.0,nonuniqueness=0.0):

        np.random.seed(seed=seed)

        sparsity = max(0,min(sparsity,1))
        nonuniqueness = max(0,min(nonuniqueness,1))
        unique_elements = int((1-sparsity)*(1-nonuniqueness)*n**2)
        repeated_elements = int((1-sparsity)*nonuniqueness*n**2)

        V = range(n)

        # Should the graph be complete with unique edge weights?

        if sparsity==0 and nonuniqueness==0:

            # If so, generate a random matrix and remove the diagonal entries.

            A = np.random.rand(n,n)
            np.fill_diagonal(A,0)

        else:

            # If not, generate a vector with the expected number of unique
            # entries, duplicate some of those entries according to a Binomial
            # distribution according to the nonuniqueness parameter, and fill
            # in the remaining entries with zeros according to the sparsity
            # parameter.  Permute the entries of the vector and reshape it into
            # the correct dimensions.

            B = np.random.rand(unique_elements)
            C = np.zeros(n**2-unique_elements,dtype=np.float)

            tally = 0
            while tally<repeated_elements:
                duplications = max(np.random.binomial(repeated_elements-tally,0.5),1)
                C[tally:tally+duplications] = B[np.random.randint(unique_elements)]
                tally += duplications

            A = np.random.permutation(np.concatenate((B,C))).reshape(n,n)
            np.fill_diagonal(A,0)

            # Is the resulting graph strongly connected?  If not, add entries
            # to the adjacency matrix at semi-random to bridge the components.

            components = strongly_connected_components(A)
            m = len(components)
            if m>1:
                for i in xrange(m):
                    for j in xrange(m):
                            if i!=j:
                                p = np.random.randint(len(components[i]))
                                q = np.random.randint(len(components[j]))
                                A[components[i][p],components[j][q]] = np.random.rand()

        return V,A

    # The next function provides the example given by Tarjan in Tarjan (1983).

    def tarjan_1983_example():

        E = [["a","b",10.0], ["b","a",12.0], ["b","c",30.0], ["d","c",6.0],  ["d","e",16.0],
             ["e","d",13.0], ["e","f",8.0],  ["f","a",26.0], ["a","g",15.0], ["g","b",35.0],
             ["c","g",45.0], ["g","c",22.0], ["d","g",14.0], ["g","e",50.0], ["f","g",20.0]]

        V,A = edges_to_matrix(E)

        return V,A,E

    # The next function tests the correctness of our implementation of the
    # hierarchical decomposition algorithm by comparing its results with those
    # from a naive implementation on various random graphs.  If the functions
    # return different trees for the same graph, then we display trees, the
    # seed for generating the adjacency matrix, and the adjacency matrix
    # itself.

    def test_correctness(m,n,increasing,sparsity=0.0,nonuniqueness=0.0):

        for i in xrange(n):

            progress("Progress: %d/%d" % (i+1,n))
            V,A = random_adjacency_matrix(m,seed=i,sparsity=sparsity,nonuniqueness=nonuniqueness)
            E = matrix_to_edges(V,A)
            S = HD_naive(E,increasing)
            T = HD(V,A,increasing)

            if S!=T:
                progress("")
                print i
                print A
                print S
                print T
                return False

        progress("")
        return True

    # The next function tests the performance of our implementation of the
    # hierarchical decomposition algorithm as well as our implementation of a
    # clustering algorithm by running then on various complete random graphs
    # of various sizes.

    def test_performance(trials,repetitions,increasing):

        import time

        dimensions = np.zeros(trials, dtype=np.int)
        HD_times = np.zeros((trials,repetitions),dtype=np.float)
        clustering_times = np.zeros((trials,repetitions),dtype=np.float)

        for trial in xrange(trials):

            n = int(10*2**trial)
            dimensions[trial] = n

            progress("Progress: %d/%d of %d/%d" % (0,repetitions,trial+1,trials))
            np.random.seed(seed=trial)
            V,A = random_adjacency_matrix(n)

            for repetition in xrange(repetitions):

                progress("Progress: %d/%d of %d/%d" % (repetition+1,repetitions,trial+1,trials))
                first_time = time.time()
                T = HD(V,A,increasing)
                second_time = time.time()
                weights,clusters = cluster(T)
                third_time = time.time()

                HD_times[trial,repetition] = second_time-first_time
                clustering_times[trial,repetition] = third_time-second_time

        progress("")

        return dimensions,np.min(HD_times,axis=1),np.min(clustering_times,axis=1)

    # The next function profiles our implementation of the hierarchical
    # decomposition algorithm on a complete random graph; search online for
    # cProfile for details.

    def test_profile(n,increasing,filename):

        import cProfile

        V,A = random_adjacency_matrix(n)
        cProfile.runctx("T = HD(V,A,"+str(increasing)+")", globals(), locals(), filename=filename)

    ###########################################################################
    #
    #   Tests
    #
    ###########################################################################

    # Excluding for the profiling code, this collection of tests should require
    # less than 1 minute of runtime and 300 MB of memory on a modern machine.
    # No additional storage space or a network connection is necessary.  We use
    # the NetworkX package only for the test cases.  We import several Fortran
    # routines to expedite several of the computations, but the code falls back
    # on functionally equivalent Python functions when the Fortran routines are
    # unavailable.

    # See Tarjan (1983) for a description of hierarchical clustering with
    # strongly connected components (SCCs).  In our implementation, we
    # represent the hierarchical decomposition (HD) tree of a strongly
    # connected graph as a dictionary.  The formation of a SCC is represented
    # as a tuple whose first entry is the condensing weight and remaining
    # entries are strongly connected vertices.

    # For example, in the example in Tarjan (1983), the addition of weight 12
    # forms SCC consisting of vertices a and b, which we represent with the
    # dictionary entries
    #
    #   (0.0, 'a') : (12.0, 'a', 'b')
    #   (0.0, 'b') : (12.0, 'a', 'b')
    #
    # The further addition of weight 35 forms a SCC consisting of vertices a,
    # b, and g, which we represent with the dictionary entries
    #
    #   (0.0, 'g') : (35.0, 'a', 'b', 'g')
    #   (12.0, 'a', 'b') : (35.0, 'a', 'b', 'g')

    # First, we run our implementation of Tarjan's HD algorithm on the example
    # given in Tarjan (1983) and output the tree, the condensing weights, and
    # the vertex clusters.  We reverse the edge order from increasing weights
    # to decreasing weights and repeat.

    print '=== Test 1 ==='
    V, A, E = tarjan_1983_example()

    T = HD(V,A,increasing=True)
    weights, clusters = cluster(T)

    print 'Weights added in increasing order:'
    print ''
    print 'Tree:'
    print prettify(T)
    print 'Condensing weights:'
    print prettify(weights)
    print 'Vertex clusters:'
    print prettify(clusters)

    T = HD(V,A,increasing=False)
    weights, clusters = cluster(T)

    print 'Weights added in decreasing order:'
    print ''
    print 'Tree:'
    print prettify(T)
    print 'Condensing weights:'
    print prettify(weights)
    print 'Vertex clusters:'
    print prettify(clusters)

    # Second, we compare tree from of our implementation of Tarjan's HD
    # algorithm with the tree from our implementation of the naive algorithm,
    # i.e., the addition of one weight at a time.  Our dictionary
    # representation of the tree is unique, so the the trees are the same when
    # the dictionaries are the same and different when different.

    print '=== Test 2 ==='
    print HD(V,A,True)==HD_naive(E,True)
    print HD(V,A,False)==HD_naive(E,False)
    print ''

    # Third, we repeat the previous test on a complete graph, a graph with
    # roughly 70% of its edges removed, and a graph with roughly 40% of its
    # edge weights repeated, many to the same weights.  Each graph has 25
    # vertices and edge weights chosen randomly from a uniform distribution.
    # Each test repeats 20 times with a different graph each time for both
    # increasing and decreasing edge weights.  The test returns `True` if the
    # trees are the same and `False` with additional information if they
    # differ; see `test_HD_correctness` for more details.

    print '=== Test 3 ==='
    print test_correctness(25,20,True)
    print test_correctness(25,20,True,sparsity=0.7)
    print test_correctness(25,20,True,nonuniqueness=0.4)
    print test_correctness(25,20,False)
    print test_correctness(25,20,False,sparsity=0.7)
    print test_correctness(25,20,False,nonuniqueness=0.4)
    print ''

    # Fourth, we generate profiling data for our implementation of Tarjan's HD
    # algorithm for a complete random graph with 2000 vertices.  Search for
    # cProfile online for details.

    print '=== Test 4 ==='
    #    test_profile(2000,True,'tarjan_HD.cprof')
    print ''

    # Fifth, we examine the performance of our implementation of Tarjan's HD
    # algorithm and the performance of our clustering routine on complete
    # random graphs of various sizes, returning the best runtime from three
    # repetitions on the same random graph.  We start with 10 vertices, double
    # to 20 vertices, and so on for a total of eight trials.  We return the
    # number of vertices and the shortest runtimes for each trial, and we
    # perform the test for both increasing and decreasing edge weights.

    print '=== Test 5 ==='
    number_of_vertices, HD_runtimes, clustering_runtimes = test_performance(8,3,True)

    print 'Number of vertices:'
    print number_of_vertices
    print 'Runtimes for increasing edge weights:'
    print HD_runtimes
    print 'Runtimes for vertex clustering:'
    print clustering_runtimes
    print ''

    number_of_vertices, HD_runtimes, clustering_runtimes = test_performance(8,3,False)

    print 'Number of vertices:'
    print number_of_vertices
    print 'Runtimes for decreasing edge weights:'
    print HD_runtimes
    print 'Runtimes for vertex clustering:'
    print clustering_runtimes
