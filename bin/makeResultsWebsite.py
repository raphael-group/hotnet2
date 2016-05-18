#!/usr/bin/python
import json
import os
import shutil
import sys
import os.path
sys.path.append(os.path.split(os.path.split(sys.argv[0])[0])[0])
import hotnet2
from hotnet2 import hnap, hnio, viz as hnviz
from hotnet2.constants import VIZ_INDEX, VIZ_SUBNETWORKS

def get_parser():
    description = 'Creates a website showing the subnetworks output by HotNet2.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-r', '--results_files', nargs='+', required=True,
                        help='Paths to results.json files output by HotNet2')
    parser.add_argument('-ef', '--edge_file', required=True,
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network.')
    parser.add_argument('-dsf', '--display_score_file',
                        help='Path to a tab seperated file contain a gene name in the first\
                              column and the display score for that gene in the second column\
                              of each line.')
    parser.add_argument('-dnf', '--display_name_file',
                        help='Path to a tab seperated file contain a gene name in the first\
                              column and the display name for that gene in the second column\
                              of each line.')
    parser.add_argument('-nn', '--network_name', default='Network',
                        help='Display name for the interaction network.')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory in which the website should be generated.')
    return parser

def load_mutation_heat(heat_file):
    # Load the heat file
    with open(heat_file, 'r') as IN:
        heat_file = json.load(IN)
        heat_parameters = heat_file['parameters']

    # Return blank mutation data if none was provided
    if not ('heat_fn' in heat_parameters and heat_parameters['heat_fn'] == 'load_mutation_heat'):
        return [], [], dict()

    # Load the mutation data
    samples = hnio.load_samples(heat_parameters['sample_file']) if heat_parameters['sample_file'] else None
    genes = hnio.load_genes(heat_parameters['gene_file']) if heat_parameters['gene_file'] else None
    snvs = hnio.load_snvs(heat_parameters['snv_file'], genes, samples) if heat_parameters['snv_file'] else []
    cnas = hnio.load_cnas(heat_parameters['cna_file'], genes, samples) if heat_parameters['cna_file'] else []

    if heat_parameters.get('sample_type_file'):
        with open(heat_parameters['sample_type_file']) as f:
            output['sampleToType'] = dict(l.rstrip().split() for l in f if not l.startswith("#") )
    else:
        if not samples:
            samples = set( m.sample for m in snvs ) | set( m.sample for m in cnas )
        sampleToType = dict( (s, "Cancer") for s in samples )

    return snvs, cnas, sampleToType

def run(args):
    subnetworks_file = '%s/viz_files/%s' % (str(hotnet2.__file__).rsplit('/', 1)[0], VIZ_SUBNETWORKS)

    # create output directory if doesn't exist; warn if it exists and is not empty
    outdir = args.output_directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if len(os.listdir(outdir)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")

    ks = set()
    output = dict(deltas=[], subnetworks=dict(), mutation_matrices=dict(), stats=dict())
    predictions = set()
    multipleHeatFiles = False
    for results_file in args.results_files:
        with open(results_file, 'r') as IN:
            results = json.load(IN)
            ccs = results['components']

        heat_file = json.load(open(results['parameters']['heat_file']))
        gene2heat = heat_file['heat']
        heat_parameters = heat_file['parameters']
        d_score = hnio.load_display_score_tsv(args.display_score_file) if args.display_score_file else None
        d_name = hnio.load_display_name_tsv(args.display_name_file) if args.display_name_file else dict()
        edges = hnio.load_ppi_edges(args.edge_file, hnio.load_index(results['parameters']['infmat_index_file']))
        delta = format(results['parameters']['delta'], 'g')
        output['deltas'].append(delta)

        output["subnetworks"][delta] = []
        predictions |= set( g for cc in ccs for g in cc )
        for cc in ccs:
            output['subnetworks'][delta].append(viz.get_component_json(cc, gene2heat, edges,
                                                                args.network_name, d_score, d_name))
        # Record the heat scores
        if 'geneToHeat' in output:
            if any( output['geneToHeat'][g] != h for g, h in gene2heat.iteritems() ) or len(gene2heat.keys()) != len(output['geneToHeat'].keys()):
                multipleHeatFiles = True
        output['geneToHeat'] = gene2heat

        # make oncoprints if heat file was generated from mutation data
        if 'heat_fn' in heat_parameters and heat_parameters['heat_fn'] == 'load_mutation_heat':
            samples = hnio.load_samples(heat_parameters['sample_file']) if heat_parameters['sample_file'] else None
            genes = hnio.load_genes(heat_parameters['gene_file']) if heat_parameters['gene_file'] else None
            snvs = hnio.load_snvs(heat_parameters['snv_file'], genes, samples) if heat_parameters['snv_file'] else []
            cnas = hnio.load_cnas(heat_parameters['cna_file'], genes, samples) if heat_parameters['cna_file'] else []
            output['geneToMutations'] = hnviz.get_component_json(ccs, snvs, cnas, d_name)

            # Get the samples and genes from the mutations directly if they weren't provided
            if not samples:
                samples = set( m.sample for m in snvs ) | set( m.sample for m in cnas )
            if not genes:
                genes = set( m.gene for m in snvs) | set( m.gene for m in cnas )

            for i, cc in enumerate(ccs):
                output['subnetworks'][delta][i]['coverage'] = viz.get_coverage(cc, snvs, cnas, samples)

            if heat_parameters.get('sample_type_file'):
                with open(heat_parameters['sample_type_file']) as f:
                    output['sampleToType'] = dict(l.rstrip().split() for l in f if not l.startswith("#") )
                    output['typeToSamples'] = dict((t, []) for t in set(output['sampleToType'].values()))
                    for s, ty in output['sampleToType'].iteritems():
                        output['typeToSamples'][ty].append( s )
            else:
                if not samples:
                    samples = set( m.sample for m in snvs ) | set( m.sample for m in cnas )
                output['sampleToType'] = dict( (s, "Cancer") for s in samples )
                output['typeToSamples'] = dict(Cancer=list(samples))

        output['stats'][delta] = results['statistics']
        ks |= set(map(int, results['statistics'].keys()))

    # Print a warning if there were multiple heat files referenced by
    # the results files
    if multipleHeatFiles:
        sys.stderr.write('Warning: results files used multiple heat files. Only the last heat file will be used to tabulate scores and make mutation matrices.\n')

    # Output to file
    output['predictions'] = sorted(predictions) # list of nodes found in any run
    output['ks'] = range(min(ks), max(ks)+1)
    with open('%s/subnetworks.json' % outdir, 'w') as out:
        json.dump(output, out, indent=4)

    shutil.copy(subnetworks_file, '%s/%s' % (outdir, VIZ_INDEX))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
