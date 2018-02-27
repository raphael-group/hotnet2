#!/usr/bin/python

# Load required modules
import sys, os, networkx as nx
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))

# Parse arguments
def get_parser():
    from hotnet2 import hnap
    description = 'Create gene-index and edge files.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-n', '--network_file', required=True,
                        help='Path to tab-separated file containing network interactions, where'\
                              'each line is of the form gene1\tgene2, indicating a network'\
                              'interaction between gene1 and gene 2.')
    parser.add_argument('-s', '--separator', required=False, default='\t',
                        help='Separator in network file; tabs by defaults')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Path to tab-separated file listing edges of the interaction network,'\
                              'where each row contains the indices of two genes that are connected'\
                              'in the network.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column'\
                              'and the name of the gene represented at that index in the second'\
                              'column of each line.')
    return parser

def run(args):
    G = nx.Graph()
    with open(args.network_file, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.rstrip('\n').split(args.separator)
                u, v = arrs[:2]
                if u!=v:
                    G.add_edge(u, v)

    G = max(nx.connected_component_subgraphs(G), key=len)

    indexToGene = dict((i+1, gene) for i, gene in enumerate(sorted(G.nodes())))
    geneToIndex = dict((gene, i) for i, gene in indexToGene.items())
    edgelist = [(geneToIndex[u], geneToIndex[v]) for u, v in G.edges()]

    with open(args.edgelist_file, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(u, v) for u, v in sorted(map(sorted, edgelist))))
    with open(args.gene_index_file, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(i, gene) for i, gene in sorted(indexToGene.items())))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
