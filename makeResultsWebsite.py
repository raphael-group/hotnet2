#!/usr/bin/python
import sys
import shutil
import os
import argparse
import json
from hotnet2 import hnio, hnap
from hotnet2.viz import *
from hotnet2.constants import *
import hotnet2

def parse_args(raw_args):
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
    parser.add_argument('-nn', '--network_name', default='Network',
                        help='Display name for the interaction network.')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory in which the website should be generated.')
    args = parser.parse_args(raw_args)

    return args

def run(args):
    index_file = '%s/viz_files/%s' % (hotnet2.__file__.rsplit('/', 1)[0], VIZ_INDEX)
    subnetworks_file = '%s/viz_files/%s' % (hotnet2.__file__.rsplit('/', 1)[0], VIZ_SUBNETWORKS)

    # create output directory if doesn't exist; warn if it exists and is not empty
    outdir = args.output_directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if len(os.listdir(outdir)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")

    deltas = list()

    for results_file in args.results_files:
        results = json.load(open(results_file))
        ccs = results['components']
        gene2heat = json.load(open(results['parameters']['heat_file']))['heat']
        d_score = hnio.load_display_score_tsv(args.display_score_file) if args.display_score_file else None
        edges = hnio.load_ppi_edges(args.edge_file)
        gene2index = dict((gene, index) for index, gene \
                          in hnio.load_index(results['parameters']['infmat_index_file']).iteritems())
        delta = results['parameters']['delta']

        deltas.append(delta)

        output = {"delta": delta, 'subnetworks': list()}
        for cc in ccs:
            output['subnetworks'].append(get_component_json(cc, gene2heat, edges, gene2index, args.network_name, d_score))

        # write output
        delta_dir = '%s/delta%s' % (outdir, delta)
        if not os.path.isdir(delta_dir):
            os.mkdir(delta_dir)
        out = open('%s/subnetworks.json' % delta_dir, 'w')
        json.dump(output, out, indent=4)
        out.close()

        shutil.copy(subnetworks_file, delta_dir)

    write_index_file(index_file, '%s/%s' % (outdir, VIZ_INDEX), deltas)

if __name__ == "__main__":
    run(parse_args(sys.argv[1:]))
