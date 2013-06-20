################################################################################
# Data loading functions

def load_index(index_file):
    arrs  = [l.split() for l in open(index_file)]
    return dict([(int(arr[0]), arr[1]) for arr in arrs])

def load_heat(heat_file):
    arrs  = [l.split() for l in open(heat_file)]
    return dict([(arr[0], float(arr[1])) for arr in arrs])

def load_gene_list(gene_list_file):
    return set([l.strip() for l in open(gene_list_file)])

################################################################################
# Data saving functions

"""
If no output file is given, the heat will be written to stdout
"""
import sys
def save_heat(heat, output_file=None):
    out = open(output_file, 'w') if output_file else sys.stdout

    for gene, score in heat.iteritems():
        out.write(gene + '\t' + str(score) + '\n')

    if output_file: out.close()