#!/usr/bin/env python

# Load required modules
import sys, os, argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directories', nargs=2, type=str, required=True)
args = parser.parse_args( sys.argv[1:] )

# Load the subnetworks per directory
def load_directory(directory):
	subnets, sig_counts = dict(), dict()
	for root, dirs, files in os.walk(directory):
		for d in dirs:
			if d.lower().startswith('delta'):
				delta = d.split('delta_')[1]

				with open('{}/{}/components.txt'.format(root, d), 'r') as IN:
					subnets[delta] = set( frozenset(l.rstrip('\n').split('\t')) for l in IN )
				with open('{}/{}/significance.txt'.format(root, d), 'r') as IN:
					sig_counts[delta] = sum( 1 for i, l in enumerate(IN) if i > 0 and float(l.rstrip('\n').split('\t')[-1]) < 0.01 )
					
	return subnets, sig_counts

subnets1, sig_counts1 = load_directory(args.directories[0])
subnets2, sig_counts2 = load_directory(args.directories[1])

print 'Deltas same?', all( d in subnets2 for d in subnets1 )
print 'Subnetwork comparison'
for d, ccs in subnets1.iteritems():
	print '\t- Same?', subnets2[d] == ccs
	print '\t- # k with P < 0.01', sig_counts1[d], sig_counts2[d]