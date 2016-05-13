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
    parser.add_argument('-nn', '--network_name', required=True,
                        help='Name of network.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output prefix.')
    parser.add_argument('-is', '--index_file_start_index', default=1, type=int,
                        help='Minimum index in the index file.')
    parser.add_argument('-b', '--beta', required=True, type=float,
                        help='Beta is the restart probability for the '\
                        'insulated heat diffusion process.')

    parser.add_argument('-q', '--Q', default=115, type=float,
                        help='Edge swap constant. The script will attempt Q*|E| edge swaps')
    parser.add_argument('-ps', '--permutation_start_index', default=1, type=int,
                        help='Index at which to start permutation file names.')
    parser.add_argument('-np', '--num_permutations', default=100, type=int,
                        help='Number of permuted networks to create.')

    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory.')

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
    pprfile = "{}/{}_ppr_{:g}.h5".format(args.output_dir, args.prefix, args.beta)
    margs = '-e %s -i %s -o %s -s %s -b %s -n %s' % (args.edgelist_file, args.gene_index_file, pprfile,
                                                     args.index_file_start_index, args.beta, args.network_name)
    ppr.run(ppr.get_parser().parse_args(margs.split()))

    # make permuted edge lists
    assert(args.num_permutations > 0)
    print "\nCreating edge lists for permuted networks"
    print "-------------------------------------------"
    if not os.path.exists(perm_dir): os.makedirs(perm_dir)
    pargs = '-q %s -s %s -e %s -p %s -o %s -n %s' % (args.Q, args.permutation_start_index, args.edgelist_file,
                                                     args.prefix, perm_dir, args.num_permutations)
    permute.run(permute.get_parser().parse_args(pargs.split()))

    # make permuted PPRs
    print "\nCreating PPR matrices for permuted networks"
    print "---------------------------------------------"
    for i in range(args.permutation_start_index, args.permutation_start_index + args.num_permutations):
        print "Working on permutation %s..." % i
        edgelist_file = '%s/%s_edgelist_%s' % (perm_dir, args.prefix, i)
        output_dir = '%s/%s' % (perm_dir, i)
        pprfile = "{}/{}_ppr_{:g}.h5".format(output_dir, args.prefix, args.beta)
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        pargs = '-e %s -i %s -o %s -s %s -b %s -n %s --exclude_network' % (edgelist_file, args.gene_index_file, pprfile,
		                                                             args.index_file_start_index, args.beta, args.network_name)
        ppr.run(ppr.get_parser().parse_args(pargs.split()))
        os.remove(edgelist_file)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
