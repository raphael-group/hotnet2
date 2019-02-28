#!/usr/bin/env python

# Load required modules
import sys, os, json
from itertools import product

# Load HotNet2 modules
from hotnet2 import run as hnrun, hnap, hnio, heat as hnheat, consensus_with_stats, viz as hnviz
from hotnet2.constants import ITERATION_REPLACEMENT_TOKEN, HN2_INFMAT_NAME

sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/scripts/'))
import createDendrogram as CD

def get_parser():
    description = "Helper script for simple runs of generalized HotNet2, including automated"\
                   "parameter selection."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-nf', '--network_files', required=True, nargs='*',
                        help='Path to HDF5 (.h5) file containing influence matrix and edge list.')
    parser.add_argument('-pnp', '--permuted_network_paths', required=True, default='',
                        help='Path to influence matrices for permuted networks, one path '\
                              'per network file. Include ' + ITERATION_REPLACEMENT_TOKEN + ' '\
                              'in the path to be replaced with the iteration number', nargs='*')
    parser.add_argument('-hf', '--heat_files', required=True, nargs='*',
                        help='Path to heat file containing gene names and scores. This can either'\
                              'be a JSON file created by generateHeat.py, in which case the file'\
                              'name must end in .json, or a tab-separated file containing a gene'\
                              'name in the first column and the heat score for that gene in the'\
                              'second column of each line.')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=2,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-d', '--deltas', nargs='*', type=float, default=[],
                        help='Delta value(s).')
    parser.add_argument('-np', '--network_permutations', type=int, default=100,
                        help='Number of permutations to be used for delta parameter selection.')
    parser.add_argument('-cp', '--consensus_permutations', type=int, default=0,
                        help='Number of permutations to be used for consensus statistical significance testing.')
    parser.add_argument('-hp', '--heat_permutations', type=int, default=100,
                        help='Number of permutations to be used for statistical significance testing.')
    parser.add_argument('-o', '--output_directory', required=True, default=None,
                        help='Output directory. Files results.json, components.txt, and'\
                              'significance.txt will be generated in subdirectories for each delta.')
    parser.add_argument('-c', '--num_cores', type=int, default=1,
                        help='Number of cores to use for running permutation tests in parallel. If'\
                              '-1, all available cores will be used.')
    parser.add_argument('-dsf', '--display_score_file',
                        help='Path to a tab-separated file containing a gene name in the first'\
                        'column and the display score for that gene in the second column of'\
                        'each line.')
    parser.add_argument('-dnf', '--display_name_file',
                        help='Path to a tab-separated file containing a gene name in the first'\
                        'column and the display name for that gene in the second column of'\
                        'each line.')
    parser.add_argument('--output_hierarchy', default=False, required=False, action='store_true',
                        help='Output the hierarchical decomposition of the HotNet2 similarity matrix.')
    parser.add_argument('--verbose', default=1, choices=range(5), type=int, required=False,
                        help='Set verbosity of output (minimum: 0, maximum: 5).')

    return parser

def run(args):
    # Load the network and heat files
    assert( len(args.network_files) == len(args.permuted_network_paths) )
    networks, graph_map = [], dict()
    for network_file, pnp in zip(args.network_files, args.permuted_network_paths):
        infmat, indexToGene, G, network_name = hnio.load_network(network_file, HN2_INFMAT_NAME)
        graph_map[network_name] = G
        networks.append( (infmat, indexToGene, G, network_name, pnp) )

    heats, json_heat_map, heat_map, mutation_map, heat_file_map = [], dict(), dict(), dict(), dict()
    for heat_file in args.heat_files:
        json_heat = os.path.splitext(heat_file.lower())[1] == '.json'
        heat, heat_name, mutations = hnio.load_heat_file(heat_file, json_heat)
        json_heat_map[heat_name] = json_heat
        heat_map[heat_name] = heat
        heat_file_map[heat_name] = heat_file
        mutation_map[heat_name] = mutations
        heats.append( (heat, heat_name) )

    # Run HotNet2 on each pair of network and heat files
    if args.verbose > 0:
        print '* Running HotNet2 in consensus mode...'

    single_runs, consensus, linkers, auto_deltas, consensus_stats = consensus_with_stats(args, networks, heats)

    # Output the single runs
    if args.verbose > 0:
        print '* Outputting results to file...'

    params = vars(args)
    result_dirs = []
    for (network_name, heat_name, run) in single_runs:
        # Set up the output directory and record for later
        output_dir = '%s/%s-%s' % (args.output_directory, network_name.lower(), heat_name.lower())
        result_dirs.append(output_dir)
        hnio.setup_output_dir(output_dir)

        # Output to file
        hnio.output_hotnet2_run(run, params, network_name, heat_map[heat_name], heat_name, heat_file_map[heat_name], json_heat_map[heat_name], output_dir)

        # create the hierarchy if necessary
        if args.output_hierarchy:
            hierarchy_out_dir = '{}/hierarchy/'.format(output_dir)
            if not os.path.isdir(hierarchy_out_dir): os.mkdir(hierarchy_out_dir)
            CD.createDendrogram( sim, index2gene.values(), hierarchy_out_dir, params, verbose=False)

    # Output the consensus
    hnio.output_consensus(consensus, linkers, auto_deltas, consensus_stats, params, args.output_directory)

    # Create the visualization(s). This has to be after the consensus procedure
    # is run because we want to default to the auto-selected deltas.
    if args.verbose > 0:
        print '* Generating and outputting visualization data...'

    d_score = hnio.load_display_score_tsv(args.display_score_file) if args.display_score_file else None
    d_name = hnio.load_display_name_tsv(args.display_name_file) if args.display_name_file else dict()
    for (network_name, heat_name, run), result_dir, auto_delta in zip(single_runs, result_dirs, auto_deltas):
        snvs, cnas, sampleToType = mutation_map[heat_name]
        G = graph_map[network_name]

        output = hnviz.generate_viz_json(run, G.edges(), network_name, heat_map[heat_name], snvs, cnas, sampleToType, d_score, d_name)

        with open('{}/viz-data.json'.format(result_dir), 'w') as OUT:
            output['params'] = dict(consensus=False, network_name=network_name, heat_name=heat_name, auto_delta=format(auto_delta, 'g'))
            json.dump( output, OUT )

    # Add the consensus visualization
    snvs, cnas, sampleToType = mutations
    consensus_ccs = [ d['core'] + d['expansion'] for d in consensus ]
    consensus_auto_delta = 0
    results = [[consensus_ccs, consensus_stats, consensus_auto_delta]]
    with open('{}/consensus/viz-data.json'.format(args.output_directory), 'w') as OUT:
        output = hnviz.generate_viz_json(results, G.edges(), network_name, heat, snvs, cnas, sampleToType, d_score, d_name)
        output['params'] = dict(consensus=True, auto_delta=format(consensus_auto_delta, 'g'))
        json.dump( output, OUT )

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
