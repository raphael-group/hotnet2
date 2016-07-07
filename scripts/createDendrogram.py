#!/usr/bin/python

# Load the required modules
import sys, os, json, pickle, os.path, scipy.io
sys.path.append(os.path.normpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
from hotnet2 import hotnet2 as hn
from hotnet2 import hnap, hnio, heat as hnheat
from hotnet2.hierarchy import HD, convertToLinkage, convertToNewick

# Parse arguments
def get_parser():
    description = 'Create a hierarchical decomposition of the HotNet2 similarity matrix.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-r', '--run_name', required=False, default='Hotnet2',
                        help='Name of run to appear in output files.')
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
    parser.add_argument('-in', '--infmat_name', required=False, default='PPR',
                        help='Name of matrix in MATLAB file.')
    parser.add_argument('-v', '--verbose', required=False, default=False, action='store_true',
                        help='Flag verbose output.')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory in which the hierarchy should be created.')
    return parser

# Create the dendrogram and output in .json and .pickle formats
def createDendrogram(sim, genespace, output_directory, params, verbose=False):
	if verbose: print '* Creating dendrogram...'

	# Perform the clustering
	if verbose: print '\t- Performing hierarchical clustering...'
	T = HD(genespace, sim, False)

	# Convert to SciPy linkage matrix (http://goo.gl/yd6z3V)
	if verbose: print '\t- Converting into SciPy linkage matrix format...'
	Z, labels = convertToLinkage(T)
	hierarchyOutput = dict(Z=Z, labels=labels, params=params)

	newickOutput = convertToNewick(T)

	# Output the linkage matrix to JSON and Newick file
	if verbose: print '\t- Outputting to JSON and Newick files...'
	with open('{}/hn2-hierarchy.json'.format(output_directory), "w") as out:
	    json.dump(hierarchyOutput, out)

	with open('{}/hn2-hierarchy.newick.txt'.format(output_directory), "w") as out:
	    out.write(newickOutput)

# Run
def run(args):
	# Load the input data
	if args.verbose: print '* Loading infmat and heat files...'
	infmat = hnio.load_infmat(args.infmat_file, args.infmat_name)
	full_index2gene = hnio.load_index(args.infmat_index_file)

	using_json_heat = os.path.splitext(args.heat_file.lower())[1] == '.json'
	if using_json_heat:
	    heat = json.load(open(args.heat_file))['heat']
	else:
	    heat = hnio.load_heat_tsv(args.heat_file)
	print "* Loaded heat scores for %s genes" % len(heat)

	# filter out genes not in the network
	heat = hnheat.filter_heat_to_network_genes(heat, set(full_index2gene.values()))

	# genes with score 0 cannot be in output components, but are eligible for heat in permutations
	heat, addtl_genes = hnheat.filter_heat(heat, None, False, 'There are ## genes with heat score 0')
	if args.verbose: print '* Creating similarity matrix...'
	sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, True)

	# Create and output the dendrogram
	createDendrogram( sim, index2gene.values(), args.output_directory, vars(args), args.verbose )

if __name__ == "__main__": run( get_parser().parse_args(sys.argv[1:]))
