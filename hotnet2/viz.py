import hnio
from collections import defaultdict

def generate_viz_json(results, edges, network_name, gene2heat, snvs, cnas, sampleToType, d_score, d_name):
    output = dict(deltas=[], subnetworks=dict(), stats=dict(), gene2heat=gene2heat)
    predictions = set()
    samples = sampleToType.keys()
    for ccs, stats, delta in results:
        delta = format(delta, 'g')
        output['stats'][delta] = stats
        output['subnetworks'][delta] = []
        for cc in ccs:
            output['subnetworks'][delta].append(get_component_json(cc, gene2heat, edges, network_name, d_score, d_name))
            predictions |= set(cc)

        if snvs or cnas:
            for i, cc in enumerate(ccs):
                output['subnetworks'][delta][i]['coverage'] = get_coverage(cc, snvs, cnas, samples)

    # Load the mutation data
    if snvs or cnas:
        output['geneToMutations'] = get_mutations_json(predictions, snvs, cnas, d_name)
        output['sampleToType'] = sampleToType

    return output

def get_nodes(cc, gene2heat, d_score, d_name):
    scores = d_score if d_score else gene2heat
    return [{'name': d_name.get(gene, gene), 'value': scores.get(gene, float('nan'))} for gene in cc]

def get_edges(cc, edges, networkName, d_name):
    edgeData = list()
    for i in range(len(cc)):
        for j in range(i+1, len(cc)):
            gene1 = min(cc[i], cc[j])
            gene2 = max(cc[i], cc[j])
            if (gene1, gene2) in edges or (gene2, gene1) in edges:
                edgeData.append({'source': d_name.get(gene1, gene1), 'target': d_name.get(gene2, gene2), 'categories': [networkName]})

    return edgeData

def get_component_json(cc, gene2heat, edges, networkName, d_score, d_name):
    nodes = get_nodes(cc, gene2heat, d_score, d_name)
    cc_edges = get_edges(cc, edges, networkName, d_name)

    return {'nodes': nodes, 'edges': cc_edges}

def get_mutations_json(genes, snvs, cnas, d_name):
    M = defaultdict(lambda: defaultdict(list))
    for mut in snvs + cnas:
        if mut.gene in genes:
            M[d_name.get(mut.gene, mut.gene)][mut.sample].append(mut.mut_type)

    return M

def get_coverage(cc, snvs, cnas, samples):
    coverage = len(set( m.sample for m in snvs + cnas if m.gene in cc ))
    return '{} ({:.2f}%)'.format(coverage, float(100*coverage)/len(samples))

def write_index_file(index_file, out_file, deltas):
    index = hnio.load_file(index_file)
    index += '<ul>\n'
    for delta in sorted(deltas):
        index += '<li><a href="delta%s/subnetworks.html">&delta; = %s</a></li>\n' % (delta, delta)
    index += '</ul>'
    hnio.write_file(out_file, index)
