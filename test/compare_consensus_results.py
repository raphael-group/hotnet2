#!/usr/bin/env python

# Load modules.
import os, sys, argparse

# Parse arguments.
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--predicted_consensus_results', type=str, required=True)
    parser.add_argument('-r', '--reference_consensus_results', type=str, required=True)
    return parser

def is_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def load_consensus_file(consensus_file):
    consensus_genes = set()
    with open(consensus_file, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                l = l.replace('[', '').replace(']', '').replace(',', '')
                genes = l.strip().split()
                for gene in genes:
                    if gene and not is_number(gene):
                        consensus_genes.add(gene)
    return consensus_genes

def run(args):
    # Load results.
    predicted_results = load_consensus_file(args.predicted_consensus_results)
    reference_results = load_consensus_file(args.reference_consensus_results)

    # Compare results.
    both_results = set.intersection(predicted_results, reference_results)
    only_predicted_results = set.difference(predicted_results, reference_results)
    only_reference_results = set.difference(reference_results, predicted_results)

    # Output results.
    print 'Summary of consensus results:'
    print '- {} genes in predicted consensus results'.format(len(predicted_results))
    print '- {} genes in reference consensus results'.format(len(reference_results))
    print '- {} genes in both predicted and reference consensus results'.format(len(both_results))
    print '- {} genes only in predicted consensus results'.format(len(only_predicted_results))
    if only_predicted_results:
        print '\t' + ', '.join(sorted(only_predicted_results))
    print '- {} genes only in reference consensus results'.format(len(only_reference_results))
    if only_reference_results:
        print '\t' + ', '.join(sorted(only_reference_results))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
