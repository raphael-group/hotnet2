import random

################################################################################
# Heat permutation

def heat_permutation_wrapper((heat_scores, eligible_genes)):
    permuted_genes = list(eligible_genes)
    random.shuffle(permuted_genes)
    permuted_genes = permuted_genes[:len(heat_scores)]

    permuted_heat = dict([(gene, heat) for gene, heat in zip(permuted_genes, heat_scores)])

    return permuted_heat

"""
heat: dict of gene name to heat score
addtl_genes: list of gene names
num_permutations: integer

Returns: list of num_permutations dicts of gene names to heat scores
"""
import multiprocessing as mp
def permute_heat(heat, num_permutations, addtl_genes=None, parallel=True):
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map

    heat_scores = heat.values()
    genes_eligible_for_heat = set(heat.keys()) | (set(addtl_genes) if addtl_genes else set())
    
    args = [(heat_scores, genes_eligible_for_heat)] * num_permutations
    permutations = map_fn(heat_permutation_wrapper, args)
    
    if parallel:
        pool.close()
        pool.join()

    return permutations

################################################################################
# Mutation permutation

def get_cna_blocks_for_sample(cnas, gene2chromo, chromo2genes):
    #sort cnas by chromosome, then by position in chromosome
    cnas.sort(key = lambda cna: (gene2chromo[cna.gene],
                                 chromo2genes[gene2chromo[cna.gene]].index(cna.gene)))
    
    chromo2blocks = {}
    curr_block = [cnas[0]]
    chromo2blocks[gene2chromo[cnas[0].gene]] = [curr_block]
    for cna in cnas[1:]:
        last_cna = curr_block[-1]
        curr_chromo = gene2chromo[last_cna.gene]
        if (gene2chromo[cna.gene] == curr_chromo and
            chromo2genes[curr_chromo].index(cna.gene) ==
                chromo2genes[curr_chromo].index(last_cna.gene) + 1 and
            cna.mut_type == last_cna.mut_type):
            curr_block.append(cna)
        else:
            curr_block = [cna]
            if gene2chromo[cna.gene] not in chromo2blocks:
                chromo2blocks[gene2chromo[cna.gene]] = []
            chromo2blocks[gene2chromo[cna.gene]] = curr_block
            
    return chromo2blocks

# def generate_random_snvs(genes_to_mutate, gene2length, bmr, gene2bmr={}):
#     """Generate a random set of SNVs and return a set of genes that are mutated as a result.
#      
#     Keyword arguments:
#     genes_to_mutate -- iterable of genes for which random SNVs should be generated
#     gene2length -- dict mapping a gene name to the gene's length
#     bmr -- default background mutation rate used for genes that do not have a gene-specific BMR
#            specified in gene2bmr
#     gene2bmr -- dict mapping a gene name to the background mutation rate for that gene. If the
#                 dict does not contain an entry for a gene in gene_to_mutate, the default BMR will
#                 be used for that 
#      
#     """
#      
#     mutated_genes = list()    
#     for gene in genes_to_mutate:
#         length = gene2length[gene]
#         if length <= 0:
#             raise RuntimeError("Non-positive length for gene %s" % gene)
#         prob = 1 - pow( 1 - bmr, length)
#         if random.random() <= prob: mutated_genes.append(gene)
#     return set(mutated_genes)

