import json
import numpy as np
import os
import shutil
import scipy.io
import sys
import hotnet2
from hotnet2 import findThreshold as ft, heat as hnheat, hnap, hnio, hotnet2 as hn, permutations as p, stats, viz
from hotnet2.constants import NUM_CCS, HEAT_JSON, JSON_OUTPUT, COMPONENTS_TSV, SIGNIFICANCE_TSV, VIZ_INDEX, VIZ_SUBNETWORKS

MIN_CC_SIZE = 3
MAX_CC_SIZE = 25
INFMAT_NAME = "Li"

def get_parser():
    description = "Helper script for simple runs of generalized HotNet2, including automated\
                   parameter selection."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-if', '--infmat_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='Path to heat file containing gene names and scores. This can either\
                              be a JSON file created by generateHeat.py, in which case the file\
                              name must end in .json, or a tab-separated file containing a gene\
                              name in the first column and the heat score for that gene in the\
                              second column of each line.')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=2,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-n', '--num_permutations', type=int, default=100,
                        help='Number of permutations that should be used for parameter selection\
                              and statistical significance testing.')
    parser.add_argument('-o', '--output_directory', default='hotnet_output',
                        help='Output directory. Files results.json, components.txt, and\
                              significance.txt will be generated in subdirectories for each delta.')
    parser.add_argument('--parallel', dest='parallel', action='store_true',
                        help='Run permutation tests in parallel.')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                        help='Run permutation tests sequentially.')
    parser.add_argument('-ef', '--edge_file',
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network. This is used to create subnetwork visualizations; if not\
                              provided, visualizations will not be made.')
    parser.add_argument('-dsf', '--display_score_file',
                        help='Path to a tab-separated file containing a gene name in the first\
                        column and the display score for that gene in the second column of\
                        each line.')
    parser.add_argument('-nn', '--network_name', default='Network',
                        help='Display name for the interaction network. (Used for subnetwork\
                              visualizations)')
    parser.set_defaults(parallel=False)
    
    return parser

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    if len(os.listdir(args.output_directory)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")
    
    infmat = scipy.io.loadmat(args.infmat_file)[INFMAT_NAME]
    full_index2gene = hnio.load_index(args.infmat_index_file)
    
    using_mutation_data = False
    using_json_heat = os.path.splitext(args.heat_file.lower())[1] == '.json'
    if using_json_heat:
        heat_data = json.load(open(args.heat_file))
        heat = heat_data['heat']
        heat_params = heat_data['parameters']
        using_mutation_data = 'heat_fn' in heat_params and heat_params['heat_fn'] == 'load_mutation_heat'
    else:
        heat = hnio.load_heat_tsv(args.heat_file)
    print "* Loading heat scores for %s genes" % len(heat)
    
    # filter out genes not in the network
    heat = hnheat.filter_heat_to_gene_set(heat, set(full_index2gene.values()), "not in the network")
    
    # genes with score 0 cannot be in output components, but are eligible for heat in permutations
    heat, addtl_genes = hnheat.filter_heat(heat, None, False, 'There are ## genes with heat score 0')

    # find delta that maximizes # CCs of size >= MIN_SIZE for each permuted data set
    deltas = ft.get_deltas_for_heat(infmat, full_index2gene, heat, addtl_genes, args.num_permutations,
                                    NUM_CCS, [MIN_CC_SIZE], True, args.parallel)

    # find the multiple of the median delta s.t. the size of the largest CC in the real data
    # is <= MAX_CC_SIZE
    medianDelta = np.median(deltas[MIN_CC_SIZE])

    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, False)

    for i in range(1, 11):
        G = hn.weighted_graph(sim, index2gene, i*medianDelta)
        max_cc_size = max([len(cc) for cc in hn.connected_components(G)])
        if max_cc_size <= MAX_CC_SIZE:
            break

    # load interaction network edges and determine location of static HTML files for visualization
    if args.edge_file:
        edges = hnio.load_ppi_edges(args.edge_file, full_index2gene)
    index_file = '%s/viz_files/%s' % (str(hotnet2.__file__).rsplit('/', 1)[0], VIZ_INDEX)
    subnetworks_file = '%s/viz_files/%s' % (str(hotnet2.__file__).rsplit('/', 1)[0], VIZ_SUBNETWORKS)

    # and run HotNet with that multiple and the next 4 multiples
    run_deltas = [i*medianDelta for i in range(i, i+5)]
    for delta in run_deltas: 
        # create output directory
        delta_out_dir = args.output_directory + "/delta_" + str(delta)
        if not os.path.isdir(delta_out_dir):
            os.mkdir(delta_out_dir)
        
        # find connected components
        G = hn.weighted_graph(sim, index2gene, delta, directed=False)
        ccs = hn.connected_components(G, args.min_cc_size)
        
        # calculate significance (using all genes with heat scores)
        print "* Performing permuted heat statistical significance..."
        heat_permutations = p.permute_heat(heat, full_index2gene.values(), args.num_permutations,
                                           addtl_genes, args.parallel)
        sizes = range(2, 11)
        print "\t- Using no. of components >= k (k \\in",
        print "[%s, %s]) as statistic" % (min(sizes), max(sizes))
        sizes2counts = stats.calculate_permuted_cc_counts(infmat, full_index2gene, heat_permutations,
                                                          delta, sizes, False, args.parallel)
        real_counts = stats.num_components_min_size(G, sizes)
        size2real_counts = dict(zip(sizes, real_counts))
        sizes2stats = stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)
    
        # sort ccs list such that genes within components are sorted alphanumerically, and components
        # are sorted first by length, then alphanumerically by name of the first gene in the component
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)
    
        # write output
        if not using_json_heat:
            heat_dict = {"heat": heat, "parameters": {"heat_file": args.heat_file}}
            heat_out = open(os.path.abspath(delta_out_dir) + "/" + HEAT_JSON, 'w')
            json.dump(heat_dict, heat_out, indent=4)
            heat_out.close()
            args.heat_file = os.path.abspath(delta_out_dir) + "/" + HEAT_JSON
        
        args.delta = delta  # include delta in parameters section of output JSON
        output_dict = {"parameters": vars(args), "sizes": hn.component_sizes(ccs),
                       "components": ccs, "statistics": sizes2stats}
        hnio.write_significance_as_tsv(os.path.abspath(delta_out_dir) + "/" + SIGNIFICANCE_TSV,
                                       sizes2stats)
        
        json_out = open(os.path.abspath(delta_out_dir) + "/" + JSON_OUTPUT, 'w')
        json.dump(output_dict, json_out, indent=4)
        json_out.close()
        
        hnio.write_components_as_tsv(os.path.abspath(delta_out_dir) + "/" + COMPONENTS_TSV, ccs)
        
        # write visualization output if edge file given
        if args.edge_file:
            viz_data = {"delta": delta, 'subnetworks': list()}
            d_score = hnio.load_display_score_tsv(args.display_score_file) if args.display_score_file else None
            for cc in ccs:
                viz_data['subnetworks'].append(viz.get_component_json(cc, heat, edges,
                                                                      args.network_name, d_score))
                
            if using_mutation_data:
                viz_data['oncoprints'] = list()
                samples = hnio.load_samples(heat_params['sample_file']) if heat_params['sample_file'] else None
                genes = hnio.load_genes(heat_params['gene_file']) if heat_params['gene_file'] else None
                snvs = hnio.load_snvs(heat_params['snv_file'], genes, samples) if heat_params['snv_file'] else []
                cnas = hnio.load_cnas(heat_params['cna_file'], genes, samples) if heat_params['cna_file'] else []
            
                for cc in ccs:
                    viz_data['oncoprints'].append(viz.get_oncoprint_json(cc, snvs, cnas))
            
            delta_viz_dir = '%s/viz/delta%s' % (args.output_directory, delta)
            if not os.path.isdir(delta_viz_dir):
                os.makedirs(delta_viz_dir)
            viz_out = open('%s/subnetworks.json' % delta_viz_dir, 'w')
            json.dump(viz_data, viz_out, indent=4)
            viz_out.close()
    
            shutil.copy(subnetworks_file, delta_viz_dir)
    
    if args.edge_file:
        viz.write_index_file(index_file, '%s/viz/%s' % (args.output_directory, VIZ_INDEX), run_deltas)
    
    
if __name__ == "__main__": 
    run(get_parser().parse_args(sys.argv[1:]))