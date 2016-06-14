#!/usr/bin/env python

# Load required modules
import sys, os, argparse, h5py, json, shutil, numpy as np

# Load required modules
parser = argparse.ArgumentParser()
parser.add_argument('-nf', '--network_files', type=str, required=True, nargs='*')
parser.add_argument('-hf', '--heat_files', type=str, required=True, nargs='*')
parser.add_argument('-o', '--output_directory', type=str, required=True)
args = parser.parse_args( sys.argv[1:] )

# Convert the network files
network_dir = args.output_directory + '/networks'
if not os.path.exists(network_dir): os.makedirs(network_dir)
for network_file in args.network_files:
	# Load the file
	f            = h5py.File(network_file, 'r')
	network      = { key:f[key].value for key in f }
	f.close()

	network_name = network['network_name']
	nodes        = network['nodes']
	edges        = network['edges']
	pnp_path     = network['permuted_networks_path']
	PPR          = np.asarray(network['PPR'])

	# Set up output directory
	this_network_dir = network_dir + '/' + network_name.split('_')[1]
	original_dir = this_network_dir + '/original'
	if not os.path.exists(original_dir): os.makedirs(original_dir)

	# Output infmat
	f        = h5py.File('{}/{}_ppr_0.4.h5'.format(original_dir, network_name), 'a')
	f['PPR'] = PPR
	f.close()
	
	# Output gene index
	nodeToIndex = dict(zip(nodes, range(1, len(nodes)+1)))
	with open('{}/{}_index_genes'.format(original_dir, network_name), 'w') as OUT:
		OUT.write('\n'.join([ '{} {}'.format(nodeToIndex[n], n) for n in nodes ]))

	# Output edge list
	with open('{}/{}_edge_list'.format(original_dir, network_name), 'w') as OUT:
		OUT.write('\n'.join(['{} {} 1'.format(nodeToIndex[u], nodeToIndex[v]) for u, v in edges ]))

	# Copy permuted networks
	permuted_dir = this_network_dir + '/permuted'
	if os.path.exists(permuted_dir): shutil.rmtree(permuted_dir)
	shutil.copytree(os.path.split(pnp_path)[0], permuted_dir)


# Copy the heat files
heat_dir = args.output_directory + '/weights'
if not os.path.exists(heat_dir): os.makedirs(heat_dir)
for heat_file in args.heat_files:
	src       = os.path.abspath(heat_file)
	file_name = os.path.split(src)[1]
	num       = file_name.split('vertex_weight_')[1].split('.')[0]
	dst_dir   =  '{}/{}'.format(heat_dir, num)
	if not os.path.exists(dst_dir): os.makedirs(dst_dir)
	shutil.copyfile(src, dst_dir + '/vertex_weight.txt')