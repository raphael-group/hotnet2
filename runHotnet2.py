# -*- coding: iso-8859-1 -*-
import hnap
import hnio
import hotnet2 as hn
import sys

def parse_args(raw_args):  
    
    description = "Runs generalized HotNet2.\
                   Note that some or all parameters can be specified via a configuration\
                   file by passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python runHotnet2.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

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
    parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    parser.add_argument('--classic', default=False, action='store_true',
                        help='Run classic (instead of directed) HotNet.')
    
    return parser.parse_args(raw_args)

def run(args):
    import scipy.io
    import json
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
    
    #infmat_index is a dict from indices to gene names
    infmat_index = hnio.load_index(args.infmat_index_file)
  
    #heat is a dict from gene names to heat scores
    heat = hnio.load_heat(args.heat_file)
  
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h, not args.classic)
    G = hn.weighted_graph(sim, gene_index, args.delta, not args.classic)
    ccs = hn.connected_components(G, args.min_cc_size)
    
    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump({"parameters": vars(args), "sizes": hn.component_sizes(ccs), "components": ccs}, output_file, indent=4)
    if (args.output_file): output_file.close()

if __name__ == "__main__": 
    run( parse_args(sys.argv[1:]) )
