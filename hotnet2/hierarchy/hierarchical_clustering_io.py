def linkage(T):
    """
    This function converts a tree from our current format to a linkage matrix Z
    and leaf node list V.
    """
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

    Y = [[a,b,base-c,d] for (a,b,c,d) in Z]
    reordered_Y, reordered_V = reorder(Y,V)

    return reordered_Y, reordered_V

def newick(T):
    """
    The next function converts a tree from our current format to the standard
    Newick format.
    """
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

def reorder(Z, V):
    """
    This function is a helper for the linkage function that reorders the leaf
    nodes to prevent crossings in the dendrogram.
    """
    # Check if the linkage matrix is constructed properly.

    for i in range(len(Z)-1):
        if Z[i][2]>Z[i+1][2]:
            raise Warning("The rows in the linkage matrix are not monotonically nondecreasing by distance; crossings may occur in the dendrogram.")
            break

    # Below, we have the following variables:

    #     m: index for reordered leaf nodes
    #     n: number of leaf nodes
    #     queue: list of inner nodes, which we append and pop as we work from right to left and from top to bottom in the dendrogram.
    #     order: dictionary for the reordering of the leaf nodes
    #     transpose_order: transpose the key, value pairs in order

    m = 0
    n = len(V)
    queue = [len(Z)-1]
    order = {}
    transpose_order = {}

    # Continue while inner nodes remain.

    while len(queue):

        # Consider the right-most inner node.

        i = queue.pop()
        subqueue = []

        # For the children of the inner node, reorder them if they are leaf nodes and append them if they are inner nodes.

        for j in reversed(range(2)):
            if Z[i][j]<n:
                order[Z[i][j]] = n-m-1
                transpose_order[n-m-1] = Z[i][j]
                m += 1
            else:
                subqueue.append(Z[i][j]-n)

        subqueue = sorted(subqueue, key=lambda i:Z[i][2])
        queue.extend(subqueue)

    # Construct the linkage matrix and list of leaf nodes for the reordered leaf nodes.

    reordered_Z = [[y for y in z] for z in Z]

    for i in range(len(Z)):
        for j in reversed(range(2)):
            if Z[i][j]<n:
                reordered_Z[i][j] = order[Z[i][j]]

    reordered_V = [V[transpose_order[i]] for i in range(n)]

    # Return the results.

    return reordered_Z, reordered_V
