#!/usr/bin/python
import os
import sys
from hotnet2 import hnap
sys.path.append('influence_matrices')
from bin import createPPRMat as ppr
from bin import permuteNetwork as permute

def get_parser():
    description = 'Create the personalized pagerank matrix and 100 permuted PPR matrices for the\
                   given network and restart probability beta.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output prefix.')
    parser.add_argument('-is', '--index_file_start_index', default=1, type=int,
                        help='Minimum index in the index file.')
    parser.add_argument('-a', '--alpha', required=True, type=float,
                        help='Page Rank dampening factor, equal to 1-beta (where beta is the\
                              restart probability for insulated heat diffusion|process).')

    parser.add_argument('-q', '--Q', default=115, type=float,
                        help='Edge swap constant. The script will attempt Q*|E| edge swaps')
    parser.add_argument('-ps', '--permutation_start_index', default=1, type=int,
                        help='Index at which to start permutation file names.')
    parser.add_argument('-n', '--num_permutations', default=100, type=int,
                        help='Number of permuted networks to create.')

    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory.')
    parser.add_argument("--matlab", default=False, action="store_true",
                        help="Create the PPR matrix using an external call "\
                             "to a MATLAB script instead of SciPy.")
    parser.add_argument("--path_to_matlab_script", default='bin/createPPRMat.m',
                    help="Path to MATLAB script if you want to use MATLAB to"\
                         " create the PPR matrix. Change this path if you are not"\
                         " running this script in default directory.")

    return parser

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
    if len(os.listdir(args.output_dir)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")

    # make real PPR
    print "\nCreating PPR matrix for real network"
    print "--------------------------------------"
    margs = '-e %s -i %s -o %s -p %s -s %s -a %s %s --path_to_matlab_script %s' % (args.edgelist_file, args.gene_index_file, args.output_dir,
                                                        args.prefix, args.index_file_start_index, args.alpha,
                                                        '--matlab' if args.matlab else '', args.path_to_matlab_script)
    ppr.run(ppr.get_parser().parse_args(margs.split()))

    # get the output edge list and index files (for the largest connected component) for permutations
    largest_cc_edgelist_file = '%s/%s_edge_list' % (args.output_dir, args.prefix)
    largest_cc_index_file = '%s/%s_index_genes' % (args.output_dir, args.prefix)

    # make permuted edge lists
    assert(args.num_permutations > 0)
    print "\nCreating edge lists for permuted networks"
    print "-------------------------------------------"
    perm_dir = '%s/permuted' % args.output_dir
    if not os.path.exists(perm_dir): os.makedirs(perm_dir)
    pargs = '-q %s -s %s -e %s -p %s -o %s -n %s' % (args.Q, args.permutation_start_index, largest_cc_edgelist_file,
                                                     args.prefix, perm_dir, args.num_permutations)
    permute.run(permute.get_parser().parse_args(pargs.split()))

    # make permuted PPRs
    print "\nCreating PPR matrices for permuted networks"
    print "---------------------------------------------"
    for i in range(args.permutation_start_index, args.permutation_start_index + args.num_permutations):
        print "Working on permutation %s..." % i
        edgelist_file = '%s/%s_edgelist_%s' % (perm_dir, args.prefix, i)
        output_dir = '%s/%s' % (perm_dir, i)
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        pargs = '-e %s -i %s -o %s -p %s -s %s -a %s %s --path_to_matlab_script %s' % (edgelist_file, largest_cc_index_file, output_dir,
                                                            args.prefix, args.index_file_start_index, args.alpha,
                                                            '--matlab' if args.matlab else '', args.path_to_matlab_script)
        ppr.run(ppr.get_parser().parse_args(pargs.split()))
        os.remove(edgelist_file)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
