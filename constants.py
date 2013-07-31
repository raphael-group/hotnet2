from collections import namedtuple

SNV = "snv"
AMP = "amp"
DEL = "del"

Mutation = namedtuple("Mutation", ["sample", "gene", "mut_type"])