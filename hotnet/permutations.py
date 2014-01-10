import random
from collections import defaultdict
from constants import *
import heat
import hnio
import multiprocessing as mp

################################################################################
# Heat permutation

def heat_permutation_wrapper((heat_scores, eligible_genes)):
    permuted_genes = list(eligible_genes)
    random.shuffle(permuted_genes)
    permuted_genes = permuted_genes[:len(heat_scores)]

    permuted_heat = dict([(gene, heat) for gene, heat in zip(permuted_genes, heat_scores)])

    return permuted_heat

def permute_heat(heat, num_permutations, addtl_genes=None, parallel=True):
    """Return a list of num_permutation dicts, each mapping gene names to permuted heat scores.
    
    Arguments:
    heat -- dict mapping a gene name to the heat score for that gene
    addtl_genes -- iterable of names of genes that do not have heat scores in the real data but
                   which may have heat scores assigned in permutations. Defaults to None.
    parallel -- whether heat permtutations should be generated in parallel. Defaults to True.
    
    """
    if parallel:
        pool = mp.Pool(25)
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

