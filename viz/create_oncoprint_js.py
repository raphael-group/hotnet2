#!/usr/bin/python
import sys
sys.path.append('../')
import json
import hnio
from collections import defaultdict

def dump_json(dictionary):
    return 'exports.data = ' + json.dumps(dictionary, sort_keys=True) + ';'

def parse_args():
    # Parse args
    import argparse
    class Args: pass
    args = Args()
    description = 'Draws an oncoprint from given subnetwork file.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input_subnetworks', required=True,
                        help='Input file.')
    parser.add_argument('-o', '--output_file', required=True,
	                help="Output file.")
    parser.add_argument('-g', '--gene_file')
    parser.add_argument('-s', '--sample_file')
    parser.add_argument('--snv_file')
    parser.add_argument('--cna_file')
    parser.parse_args(namespace=args)

    return args

# Oncoprint visualization data
def subgraph_mutations_patterns(subgraph):
    '''Return a map from each sample to whether s/he is mutated exclusively
    in the given subgraph.'''
    h = {}
    for g, sample_data in subgraph:
        for sample_id, mutations in sample_data:

            if sample_id in h.keys(): h[sample_id] = False
            else: h[sample_id] = True
    return h

def tys_to_int(tys):
    '''Converts a list of mutation types (i.e. snv/amp/del) into an integer
    that is easily sorted.'''
    snv, delete, amp = "snv" in tys, "del" in tys, "amp" in tys
    if (snv and amp) or (snv and delete): return 0
    elif snv: return 1
    else: return 2

def sort_samples(subgraph, samples, samples2exclusivity):
    '''Return a sorted list of ALL samples using the following priority:
    - exclusve mutations first
    - by cancer type (alphabetical)
    - mutation type (snvs first (including multiple types), then cnas)
    - by sample id (alphabetical)
    '''        
    sorted_samples = []
    for g, mut_data in subgraph:
        mut_ids = set([sid for sid, ty in mut_data])
        sid2mut_type = dict([(sid, tys_to_int(ty)) for sid, ty in  mut_data])
        mut_samples = filter(lambda (sid, ty): sid in mut_ids, samples)
        mut_samples.sort(key=lambda (sid, ty):
           (int(not samples2exclusivity[sid]), ty, sid2mut_type[sid], sid))
        sorted_samples += [(sid, ty) for sid, ty in mut_samples]
        for sample, cancer_type in mut_samples:
            samples.remove((sample, cancer_type))
    return sorted_samples # + samples # (don't include samples with no mutations)

def sort_by_frequency(subgraph):
    '''Sort genes by the number of samples they are mutated in'''
    return sorted(subgraph, key=lambda (g, ids): len(ids), reverse=True)

def clean_samples(samples):
    '''Remove unncessary IDs (e.g. TCGA-) from a list of samples.'''
    return map(lambda s: s.split('TCGA-')[1], samples)

def create_oncoprint(subnetwork, gene2heat, samples_w_types):
    subnetwork_w_heat = []
    for g in subnetwork:
        subnetwork_w_heat.append( [ g, list(gene2heat[g].items()) ] )
        
    # subnetwork_w_heat = sort_by_frequency(subnetwork_w_heat)
    samples2exclusivity = subgraph_mutations_patterns(subnetwork_w_heat)
    sorted_samples = sort_samples(subnetwork_w_heat, list(samples_w_types),
                                  samples2exclusivity)
    
    mutation_matrix = []
    genes_with_cov  = []
    for g, mut_samples in subnetwork_w_heat:
        genes_with_cov.append( "%s (%s)" % (g, len(mut_samples)) )
        pattern = []
        mut_sample_ids = [sid for sid, muts in mut_samples]
        for sample, ty in sorted_samples:
            mut_data = { "snv" : False, "del" : False, "amp" : False,
                         "cancer_type" : ty, "type" : "N/A" }
            if sample in mut_sample_ids:
                for mut_ty in mut_samples[mut_sample_ids.index(sample)][1]:
                    mut_type = "ex" if samples2exclusivity[sample] else "co"
                    mut_data["type"] = mut_type
                    mut_data[mut_ty] = True
                    if mut_ty == "amp": mutation_name = g + "(A)"
                    elif mut_ty == "del": mutation_name = g + "(D)"
                    elif mut_ty == "snv": mutation_name = g
                    else: print mut_type
                    
            pattern.append(mut_data)
        mutation_matrix.append({ "gene" : g, "pattern" : pattern })

    clean_sample_names = clean_samples([s for s, ty in sorted_samples])
    return dict(mutation_matrix=mutation_matrix,
                genes=genes_with_cov,
                samples=clean_sample_names)

def run(args):    
    # Load mutation data
    print "* Loading mutation data..."
    
    if not args.snv_file and not args.cna_file:
        raise ValueError("")
    
    samples = hnio.load_samples(args.sample_file) if args.sample_file else None
    genes = hnio.load_genes(args.gene_file) if args.gene_file else None
    snvs = hnio.load_snvs(args.snv_file, genes, samples) if args.snv_file else []
    cnas = hnio.load_cnas(args.cna_file, genes, samples) if args.cna_file else []
    mutations = snvs + cnas
    
    gene2heat = defaultdict(dict)
    for mut in mutations:
        gene2heat[mut.gene][mut.sample] = [mut.mut_type]
    
    if not samples:
        samples = set([mut.sample for mut in mutations])
    samples_w_types = [(sample, "Unknown") for sample in samples]
 
    # Load subnetwork
    print "* Loading subnetworks..."
    arrs = [l.rstrip().split("\t") for l in open(args.input_subnetworks) if not l.startswith("#")]
    name2subnetwork = dict([("Subnetwork%s" % i, arrs[i]) for i in range(len(arrs))])
     
    # Draw oncoprint
    print "* Create oncoprint json files..."
    oncoprints = dict()
    for name, subnetwork in name2subnetwork.items():
        #if name != "RTK": continue
        coverage = len(set([ p for g in subnetwork for p in gene2heat[g].keys() ]))
        print "\t", name, "(%s)" % coverage
        subnetwork.sort(key=lambda g: len(gene2heat[g]), reverse=True)
        oncoprints[name] = create_oncoprint( subnetwork, gene2heat, samples_w_types )
         
    # Output JSON
    open(args.output_file, 'w').write( dump_json( oncoprints ))

if __name__ == "__main__": run( parse_args() )
