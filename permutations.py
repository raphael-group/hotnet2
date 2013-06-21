################################################################################
# Permutation functions

from random import shuffle
def heat_permutation_wrapper((heat_scores, eligible_genes)):
    permuted_genes = list(eligible_genes)
    shuffle(permuted_genes)
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



