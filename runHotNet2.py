import sys
from hotnet2 import run as hnrun, hnap
from hotnet2.constants import ITERATION_REPLACEMENT_TOKEN

MAX_CC_SIZES = [5, 10, 15, 20]
INFMAT_NAME = "PPR"

def get_parser():
    description = "Helper script for simple runs of generalized HotNet2, including automated\
                   parameter selection."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to HDF5 (.h5) file containing influence matrix. '\
                        'NumPy (.npy) and MATLAB (.mat) files also supported.')
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
    parser.add_argument('-pnp', '--permuted_networks_path', required=False, default='',
                        help='Path to influence matrices for permuted networks. Include ' +\
                              ITERATION_REPLACEMENT_TOKEN + ' in the path to be replaced with the\
                              iteration number')
    parser.add_argument('-d', '--deltas', nargs='*', type=float, default=[],
                        help='Delta value(s).')
    parser.add_argument('-dp', '--delta_permutations', type=int, default=100,
                        help='Number of permutations to be used for delta parameter selection.')
    parser.add_argument('-sp', '--significance_permutations', type=int, default=100,
                        help='Number of permutations to be used for statistical significance testing.')
    parser.add_argument('-o', '--output_directory', default='hotnet_output',
                        help='Output directory. Files results.json, components.txt, and\
                              significance.txt will be generated in subdirectories for each delta.')
    parser.add_argument('-c', '--num_cores', type=int, default=1,
                        help='Number of cores to use for running permutation tests in parallel. If\
                              -1, all available cores will be used.')
    parser.add_argument('-ef', '--edge_file',
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network. This is used to create subnetwork visualizations; if not\
                              provided, visualizations will not be made.')
    parser.add_argument('-dsf', '--display_score_file',
                        help='Path to a tab-separated file containing a gene name in the first\
                        column and the display score for that gene in the second column of\
                        each line.')
    parser.add_argument('-dnf', '--display_name_file',
                        help='Path to a tab-separated file containing a gene name in the first\
                        column and the display name for that gene in the second column of\
                        each line.')
    parser.add_argument('-nn', '--network_name', default='Network',
                        help='Display name for the interaction network. (Used for subnetwork\
                              visualizations)')
    parser.add_argument('--output_hierarchy', default=False, required=False, action='store_true',
                        help='Output the hierarchical decomposition of the HotNet2 similarity matrix.')

    return parser

def run(args):
    extra_delta_args = [args.permuted_networks_path, INFMAT_NAME, MAX_CC_SIZES]
    hnrun.run_helper(args, INFMAT_NAME, hnrun.get_deltas_hotnet2, extra_delta_args)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
