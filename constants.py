from collections import namedtuple

#mutation types
SNV = "snv"
AMP = "amp"
DEL = "del"
INACTIVE_SNV = "inactive_snv"

#test statistics
MAX_CC_SIZE = 'max_cc_size'
NUM_CCS = 'num_ccs'

#output file names
JSON_OUTPUT = "results.json"
COMPONENTS_TSV = "components.txt"
SIGNIFICANCE_TSV = "significance.txt"

Mutation = namedtuple("Mutation", ["sample", "gene", "mut_type"])

ITERATION_REPLACEMENT_TOKEN = '##NUM##'
