#!/usr/bin/python
import os
import sys
from hotnet2 import hnap
sys.path.append('influence_matrices')
import createPPRMat as ppr
import permuteNetwork as permute

def parse_args(raw_args):
    description = 'Create the personalized pagerank matrix and 100 permuted PPR matrices for the\
                   given network and restart probability beta.'
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-e', '--edgelist_file', required=True,
                        help='Location of edgelist file.')
    parser.add_argument('-i', '--gene_index_file', required=True,
                        help='Location of gene-index file.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output prefix.')
    parser.add_argument('-is', '--index_file_start_index', default=1, type=int,
                        help="Minimum index in the inde file.")
    parser.add_argument('-a', '--alpha', required=True, type=float,
                        help="Page Rank dampening factor.")

    parser.add_argument('-q', '--Q', default=115, type=float,
                        help='Edge swap constant.')
    parser.add_argument('-ps', '--permutation_start_index', default=1, type=int,
                        help='Index at which to start permutation file names.')
    parser.add_argument('-n', '--num_permutations', default=100, type=int,
                        help='Number of permuted networks to create.')

    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory.')
    parser.add_argument("--matlab", default=False, action="store_true",
                        help="Create the PPR matrix using an external call "\
                             "to a MATLAB script instead of Scipy.")

    return parser.parse_args(raw_args)

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
    if len(os.listdir(args.output_dir)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. "
              "(Ctrl-c to cancel).")

    # make real PPR
    print "\nCreating PPR matrix for real network"
    print "--------------------------------------"
    margs = '-e %s -i %s -o %s -p %s -s %s -a %s %s' % (args.edgelist_file, args.gene_index_file, args.output_dir,
                                                        args.prefix, args.index_file_start_index, args.alpha,
                                                        '--matlab' if args.matlab else '')
    ppr.run(ppr.parse_args(margs.split()))

    # get the output edge list and index files (for the largest connected component) for permutations
    largest_cc_edgelist_file = '%s/%s_edge_list' % (args.output_dir, args.prefix)
    largest_cc_index_file = '%s/%s_index_genes' % (args.output_dir, args.prefix)

    # make permuted edge lists
    print "\nCreating edge lists for permuted networks"
    print "-------------------------------------------"
    perm_dir = '%s/permuted' % args.output_dir
    if not os.path.exists(perm_dir): os.makedirs(perm_dir)
    pargs = '-q %s -s %s -e %s -p %s -o %s -n %s' % (args.Q, args.permutation_start_index, largest_cc_edgelist_file,
                                                     args.prefix, perm_dir, args.num_permutations)
    permute.run(permute.parse_args(pargs.split()))

    # make permuted PPRs
    print "\nCreating PPR matrices for permuted networks"
    print "---------------------------------------------"
    for i in range(args.permutation_start_index, args.permutation_start_index + args.num_permutations):
        print "Working on permutation %s..." % i
        edgelist_file = '%s/%s_edgelist_%s' % (perm_dir, args.prefix, i)
        output_dir = '%s/%s' % (perm_dir, i)
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        pargs = '-e %s -i %s -o %s -p %s -s %s -a %s %s' % (edgelist_file, largest_cc_index_file, output_dir,
                                                            args.prefix, args.index_file_start_index, args.alpha,
                                                            '--matlab' if args.matlab else '')
        ppr.run(ppr.parse_args(pargs.split()))
        os.remove(edgelist_file)

if __name__ == "__main__":
    run(parse_args(sys.argv[1:]))
