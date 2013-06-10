# -*- coding: iso-8859-1 -*-
from hotnet2 import *

#random note: /data/compbio/datasets/HeatKernels/IntHint/binary+complex+hi2012/inthint_inf_0.10.mat contains 2 9859x9849 matrices, L and Li
#to save memory loading, probably best to have this only contain Li

def parse_args(raw_args):  
    import argparse
    
    class BetterFileArgParser(argparse.ArgumentParser):
        def convert_arg_line_to_args(self, arg_line):
            for arg in arg_line.split():
                yield arg
            
    description = "Runs generalized HotNet2.\
                   Note that some or all parameters can be specified via a configuration\
                   file by passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python runHotnet2.py @testConf.txt --runname TestRun'."
    parser = BetterFileArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-mn', '--infmat_name', required=True, default='Li',
                        help='Variable name of the influence matrix in the .mat file')
    parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                        help='Gene-index file for the influence matrix.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='Heat score file')
    parser.add_argument('-d', '--delta', type=float, required=True,
                        help='Weight threshold for edge removal')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=3,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    
    return parser.parse_args(raw_args)

def load_index( index_file ):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_heat( heat_file ):
	arrs  = [l.split() for l in open(heat_file)]
	return dict([(arr[0], float(arr[1])) for arr in arrs])

def run(args):
	import scipy.io
	infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
	
	#infmat_index is a dict from indices to gene names
	infmat_index = load_index(args.infmat_index_file)
  
	#heat is a dict from gene names to heat scores
	heat = load_heat(args.heat_file)
  
	M, gene_index, inf_score = induce_infmat(infmat, infmat_index, sorted(heat.keys()))
	h = heat_vec(heat, gene_index)
	sim, sim_score = similarity_matrix(M, h, gene_index, not args.classic)
	G = weighted_graph(sim, gene_index, args.delta)
	
	print "* Sizes:"
	print component_sizes(G, args.min_cc_size)
	print "* Components:"
	print write_components(G, args.min_cc_size)

if __name__ == "__main__": 
	from sys import argv
	run( parse_args(argv[1:]) )