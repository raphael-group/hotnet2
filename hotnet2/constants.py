from collections import namedtuple

# Names
HOTNET2 = "hotnet2"
HOTNET = "hotnet"

# mutation types
SNV = "snv"
AMP = "amp"
DEL = "del"
INACTIVE_SNV = "inactive_snv"
FUSION = "fus"

# test statistics
MAX_CC_SIZE = 'max_cc_size'
NUM_CCS = 'num_ccs'

# output file names
HEAT_JSON = "heat.json"
JSON_OUTPUT = "results.json"
COMPONENTS_TSV = "components.txt"
SIGNIFICANCE_TSV = "significance.txt"

VIZ_INDEX = 'index.html'
VIZ_SUBNETWORKS = 'subnetworks.html'
HIERARCHY_WEB_FILE = 'hierarchy.html'

Mutation = namedtuple("Mutation", ["sample", "gene", "mut_type", "valid"])
Mutation.__new__.__defaults__ = (True,) # default valid field to True
Fusion = namedtuple("Fusion", ["sample", "genes"])

ITERATION_REPLACEMENT_TOKEN = '##NUM##'

# HotNet2 params
HN2_MAX_CC_SIZES = [5, 10, 15, 20]
HN2_INFMAT_NAME = "PPR"
HN2_STATS_SIZES = range(2, 11)
