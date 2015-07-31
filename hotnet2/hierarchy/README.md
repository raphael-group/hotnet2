Overview
-----------
Tarjan (1982, 1983) describes a hierarchical clustering algorithm that uses strongly connected components.  Below, we describe how to use Tarjan's algorithm to select the `delta` parameter in HotNet2.

Installation
------------
Run the included `setup_fortran.py` script or, (almost) equivalently, run the following command in terminal:

    f2py -c fortran_routines.f95 -m fortran_routines

This process requires a Fortran compiler and NumPy, which contains the F2PY program.  Our implementation uses Fortran in several places where large performance gains are possible.  However, if a Fortran compiler is unavailable or unsuccessful, our implementation will automatically revert to functionally equivalent (but slower) Python functions.

Usage
-----
For illustration, consider the example given in Tarjan (1983).  The graph in this example can be given by the vertices `V` and weighted adjacency matrix `A`.

`V`:

    ['a', 'b', 'c', 'd', 'e', 'f', 'g']

`A`:

    [[  0.  10.   0.   0.   0.   0.  15.]
     [ 12.   0.  30.   0.   0.   0.   0.]
     [  0.   0.   0.   0.   0.   0.  45.]
     [  0.   0.   6.   0.  16.   0.  14.]
     [  0.   0.   0.  13.   0.   8.   0.]
     [ 26.   0.   0.   0.   0.   0.  20.]
     [  0.  35.  22.   0.  50.   0.   0.]]

For example, the edge from vertex `a` to vertex `b` has weight `10.`.

For our implementation, the entries of the list `V` may be any hashable type, e.g., integers and/or strings.  The entries of the NumPy array `A` must be double-precision floating-point numbers.  The entries of `A` need not be unique.  An exception is raised when the graph described by `A` is not strongly connected.

To generate the hierarchical decomposition tree, import `hierarchical_clustering.py` and run

    T = HD(V,A,increasing)

where `increasing = True` adds edge weights in increasing order and `increasing = False` adds edge weights in decreasing order; the latter is the case for HotNet2.  The `increasing` argument is optional, and `increasing = False` is its default value when omitted.  For this example, set `increasing=True`.  The output is a tree represented as a dictionary `T` whose keys and values are tuples containing edge weights and leaf nodes.

`T`:

    (0.0, 'a') : (12.0, 'a', 'b')
    (0.0, 'b') : (12.0, 'a', 'b')
    (0.0, 'c') : (45.0, 'a', 'b', 'c', 'g')
    (0.0, 'd') : (16.0, 'd', 'e')
    (0.0, 'e') : (16.0, 'd', 'e')
    (0.0, 'f') : (50.0, 'a', 'b', 'c', 'd', 'e', 'f', 'g')
    (0.0, 'g') : (35.0, 'a', 'b', 'g')
    (12.0, 'a', 'b') : (35.0, 'a', 'b', 'g')
    (16.0, 'd', 'e') : (50.0, 'a', 'b', 'c', 'd', 'e', 'f', 'g')
    (35.0, 'a', 'b', 'g') : (45.0, 'a', 'b', 'c', 'g')
    (45.0, 'a', 'b', 'c', 'g') : (50.0, 'a', 'b', 'c', 'd', 'e', 'f', 'g')

For example, vertices `a` and `b` condense into a single component after the addition of weight `12.0`.  The later addition of `35.0` combines the component containing both `a` and `b` with the component containing only `g` into a single component containing `a`, `b`, and `g`.  The eventual addition of `50.0` completes the tree, combining `a`, `b`, `c`, `d`, `e`, `f`, and `g` into a single component.

To convert this specialized tree representation to the more standard Newick format, run

    newick_representation = newick(T)

The output is a string representing the tree.

`newick_representation`:

    (f:50.0,(d:16.0,e:16.0):34.0,(c:45.0,(g:35.0,(a:12.0,b:12.0):23.0):10.0):5.0);

To find the clusters formed from the hierarchical decomposition, run

    weights,clusters = cluster(T)

The output is a collection of condensing weights and a collection of vertex clusters.

`weights`:

     [0.0, 12.0, 16.0, 35.0, 45.0, 50.0]

`clusters`:

     [['a'], ['b'], ['c'], ['d'], ['e'], ['f'], ['g']]
     [['a', 'b'], ['c'], ['d'], ['e'], ['f'], ['g']]
     [['a', 'b'], ['c'], ['d', 'e'], ['f'], ['g']]
     [['a', 'b', 'g'], ['c'], ['d', 'e'], ['f']]
     [['a', 'b', 'c', 'g'], ['d', 'e'], ['f']]
     [['a', 'b', 'c', 'd', 'e', 'f', 'g']]

For example, each vertex initially forms its own strongly connected component.  Later, after the addition of `12.0`, the components with `a` and `b` merge to form one component with both `a` and `b`.  After the eventual addition of `50.0`, the vertices `a`, `b`, `c`, `d`, `e`, `f`, and `g` form a single connected component.

Minimal working example
-----------------------
Consider the following minimal working example.  Note that `increasing` has been omitted from `HD`, so it is `False` by default.

** Input **

    from hierarchical_clustering import *
    V = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    A = np.array([[  0.,  10.,   0.,   0.,   0.,   0.,  15.],
                  [ 12.,   0.,  30.,   0.,   0.,   0.,   0.],
                  [  0.,   0.,   0.,   0.,   0.,   0.,  45.],
                  [  0.,   0.,   6.,   0.,  16.,   0.,  14.],
                  [  0.,   0.,   0.,  13.,   0.,   8.,   0.],
                  [ 26.,   0.,   0.,   0.,   0.,   0.,  20.],
                  [  0.,  35.,  22.,   0.,  50.,   0.,   0.]])

    T = HD(V,A)
    weights,clusters = cluster(T)
    newick_representation = newick(T)

    print 'Tree in Newick format:'
        print '    ',newick_representation
    print 'Tree in our specialized format:'
    for v in T:
        print '    ',v,':',T[v]
    print 'Condensing weights:'
    print '    ',weights
    print 'Clusters:'
    for c in clusters:
        print '    ',c

** Output **

    Tree in Newick format:
        (f:42.0,(a:38.0,(d:37.0,e:37.0,(b:20.0,c:20.0,g:20.0):17.0):1.0):4.0);

    Tree in our specialized format:
        (50.0, 'a') : (12.0, 'a', 'b', 'c', 'd', 'e', 'g')
        (50.0, 'b') : (30.0, 'b', 'c', 'g')
        (50.0, 'c') : (30.0, 'b', 'c', 'g')
        (13.0, 'b', 'c', 'd', 'e', 'g') : (12.0, 'a', 'b', 'c', 'd', 'e', 'g')
        (50.0, 'd') : (13.0, 'b', 'c', 'd', 'e', 'g')
        (50.0, 'e') : (13.0, 'b', 'c', 'd', 'e', 'g')
        (50.0, 'f') : (8.0, 'a', 'b', 'c', 'd', 'e', 'f', 'g')
        (30.0, 'b', 'c', 'g') : (13.0, 'b', 'c', 'd', 'e', 'g')
        (50.0, 'g') : (30.0, 'b', 'c', 'g')
        (12.0, 'a', 'b', 'c', 'd', 'e', 'g') : (8.0, 'a', 'b', 'c', 'd', 'e', 'f', 'g')

    Condensing weights:
        0 : 30.0
        1 : 13.0
        2 : 12.0
        3 : 8.0
        4 : 0.0

    Vertex clusters:
        0 : [['a'], ['b'], ['c'], ['d'], ['e'], ['f'], ['g']]
        1 : [['a'], ['b', 'c', 'g'], ['d'], ['e'], ['f']]
        2 : [['a'], ['b', 'c', 'd', 'e', 'g'], ['f']]
        3 : [['a', 'b', 'c', 'd', 'e', 'g'], ['f']]
        4 : [['a', 'b', 'c', 'd', 'e', 'f', 'g']]

References
----------
* R.E. Tarjan. A Hierarchical Clustering Algorithm Using Strong Components. *Information Processing Letters*. (1982) 14(1):26-29.
* R.E. Tarjan. An Improved Algorithm for Hierarchical Clustering Using Strong Components. *Information Processing Letters*. (1983) 17(1):37-41.
