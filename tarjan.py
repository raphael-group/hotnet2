from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import math
import random
import string

RWL = 4
def randomword(length):
    return ''.join(random.choice(string.lowercase) for _ in range(length))

def e(node1, node2, weight):
    return (node1, node2, {"weight": weight})

def get_index(edges, weight):
    index = -1
    for i in range(len(edges)):
        if edges[i][2]['weight'] <= weight:
            index = i
    print 'index', index
    return index

def condensation(G, sccs):
    edges = sorted(G.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)
    mapping = {}
    names = defaultdict(str)
    for i, component in enumerate(sccs):
        for n in component:
            mapping[n] = i
            names[i] += n
    C = nx.DiGraph()
    C.add_nodes_from(names.values())
    for edge in edges:
        if mapping[edge[0]] != mapping[edge[1]]:
            C.add_edge(names[mapping[edge[0]]], names[mapping[edge[1]]], {'weight': edge[2]['weight']})
    return C

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def get_root(G):
    root = -1
    for node in G.nodes():
        if type(node) == int or is_int(node.split('_')[0]):
            root = node
    if root == -1:
        raise ValueError()
    return root

def tarjan(G, edges, i):
    if len(nx.strongly_connected_components(G)) != 1:
        raise ValueError("G must be strongly connected")
    
    print '\n', G.nodes(), ' - ', edges, ' - ', i
    if len(edges) - (i+1) == 1 or len(edges) - (i+1) == 0:
        print "in base case"
        root = edges[-1][2]['weight']
        #root = str(edges[-1][2]['weight']) + '_' + randomword(RWL)
        T = nx.DiGraph()
        T.add_edges_from([(root, node) for node in G.nodes()])
        return T
    else:
        j = int(math.ceil(((i+1) + len(edges)) / 2.0))
        print j
        G_j = nx.DiGraph()
        G_j.add_nodes_from(G.nodes())
        G_j.add_edges_from(edges[:j])
        print "Gj nodes", G_j.nodes()
        print "Gj edges", G_j.edges(data=True)
        scc_graphs = nx.strongly_connected_component_subgraphs(G_j)
        if len(scc_graphs) == 1:
            return tarjan(G_j, edges[:j], i)
        else:
            subgraph_trees = {}
            for scc_graph in scc_graphs:
                if len(scc_graph) > 1:
                    scc_edges = sorted(scc_graph.edges(data=True), key=lambda x: x[2]['weight'])
                    scc_index = get_index(scc_edges, edges[i][2]['weight']) if i > -1 else -1
                    subgraph_trees[''.join(scc_graph.nodes())] = tarjan(scc_graph, scc_edges, scc_index)
        C = condensation(G, [scc_graph.nodes() for scc_graph in scc_graphs])
        print C.nodes()
        print C.edges(data=True)
        c_edges = sorted(C.edges(data=True), key=lambda x: x[2]['weight'])
        c_index = get_index(c_edges, edges[j-1][2]['weight'])
        T_c = tarjan(C, c_edges, c_index)
        root = get_root(T_c)
        print root
        to_return = nx.DiGraph()
        print "keys:", subgraph_trees.keys()
        for edge in T_c.edges():
            if edge[1] in subgraph_trees:
                subgraph_tree = subgraph_trees[edge[1]]
                to_return.add_edge(edge[0], get_root(subgraph_tree))
                to_return.add_edges_from(subgraph_tree.edges())
            else:
                to_return.add_edge(edge[0], edge[1])
        print "to_return:", to_return.nodes(), to_return.edges()
        return to_return

def draw(G):
    nx.draw(G)
    plt.show()

# testEdges = [e('a', 'b', 1), e('b', 'a', 2)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 2), e('b', 'c', 3), e('c', 'a', 4)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1), e('b', 'c', 2), e('c', 'b', 2)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1), e('c', 'd', 2), e('d', 'c', 2),
#             e('a', 'c', 3), e('c', 'a', 3)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1), e('c', 'd', 1), e('d', 'c', 1),
#             e('a', 'c', 3), e('c', 'a', 3)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1), e('c', 'd', 2), e('d', 'c', 2),
#              e('e', 'f', 3), e('f', 'e', 3), e('a', 'c', 4), e('c', 'a', 4),
#              e('e', 'c', 5), e('c', 'e', 5)]

# testEdges = [e('a', 'b', 1), e('b', 'a', 1), e('c', 'd', 1), e('d', 'c', 1),
#              e('e', 'f', 1), e('f', 'e', 1), e('a', 'c', 3), e('c', 'a', 3),
#              e('e', 'd', 3), e('d', 'e', 3)]

#########################

# testEdges = [e('a', 'b', 10),  e('a', 'g', 15), e('b', 'a', 12),  e('g', 'b', 35), e('b', 'c', 30),
#              e('c', 'g', 45)]

testEdges = [e('a', 'b', 10),  e('a', 'g', 15), e('b', 'a', 12), e('b', 'c', 30), e('c', 'g', 45),
             e('d', 'c', 6), e('d', 'e', 16), e('d', 'g', 14), e('e', 'd', 13), e('e', 'f', 8),
             e('f', 'a', 26), e('f', 'g', 20), e('g', 'b', 35), e('g', 'c', 22), e('g', 'e', 50)]

#########################

# testEdges = [e('a', 'b', 10), e('b', 'a', 10), e('b', 'c', 5), e('c', 'a', 5)]

# testEdges = [e('a', 'b', -10), e('b', 'a', -10), e('b', 'c', -5), e('c', 'a', -5)]

testG = nx.DiGraph()
testG.add_edges_from(testEdges)

# draw(testG)

T = tarjan(testG, sorted(testEdges, key=lambda x: x[2]['weight']), -1)
draw(T)