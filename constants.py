from collections import namedtuple

#mutation types
SNV = "snv"
AMP = "amp"
DEL = "del"

#test statistics
MAX_CC_SIZE = 'max_cc_size'
NUM_CCS = 'num_ccs'

Mutation = namedtuple("Mutation", ["sample", "gene", "mut_type"])