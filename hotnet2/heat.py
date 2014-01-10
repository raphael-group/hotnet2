from math import log10
from collections import defaultdict
import scipy
from constants import *

def filter_heat(heat, min_score):
    """Returns (1) a dict mapping gene names to heat scores that contains only entries with heat
    scores greater than min_score; (2) a list of the genes excluded because their heat scores were
    less than min_score; and (3) the value of min_score, which will not be None even if
    it was given as None.
    
    Arguments:
    heat -- dict mapping gene names to heat scores
    min_score -- minimum heat score a gene must have to be included in the returned dict. If None,
                 the minimum score will be calculated as the minimum of the non-zero heat scores
                 included in the input dict 
    """
    if min_score is None:
        min_score = min([score for score in heat.values() if score > 0])
    filtered_heat = dict([(gene, score) for gene, score in heat.iteritems() if score >= min_score])
    return filtered_heat, [gene for gene in heat if gene not in filtered_heat], min_score

def num_snvs(mutations):
    """Return the number of SNVs in the given iterable of Mutations (those with mut_type == SNV).
    
    Arguments:
    mutations -- iterable of Mutation tuples 
    
    """
    return len([mut for mut in mutations if mut.mut_type == SNV])

def num_cnas(mutations):
    """Return the number of CNAs in the given iterable of Mutations (those with mut_type == AMP or DEL).
    
    Arguments:
    mutations -- iterable of Mutation tuples 
    
    """
    return len([mut for mut in mutations if mut.mut_type == AMP or mut.mut_type == DEL])

def filter_cnas(cnas, filter_thresh):
    """Return a list of Mutation tuples containing only CNAs in genes in which the proportion of
    CNAs in the gene across samples is at least filter_thresh. CNAs in genes whose CNAs pass the
    threshold that are of the opposite type of the dominant type in the gene are also excluded.
    
    Arguments:
    cnas -- a list of Mutation tuples representing CNAs. All are assumed to be of mut_type AMP or DEL
    filter_thresh -- the proportion of CNAs in a gene that must be of the same type
    
    """
    if filter_thresh <= .5:
        raise ValueError("filter_thresh must be greater than .5")
    
    filtered_cnas = list(cnas)
    genes2cnas = defaultdict(list)
    for cna in filtered_cnas:
        genes2cnas[cna.gene].append(cna)
    
    for gene_cnas in genes2cnas.itervalues():
        amp_count = float(len([cna for cna in gene_cnas if cna.mut_type == AMP]))
        del_count = float(len([cna for cna in gene_cnas if cna.mut_type == DEL]))
        if (amp_count / (amp_count + del_count)) >= filter_thresh:
            remove_opposite_cnas(filtered_cnas, gene_cnas, AMP)
        elif (del_count / (amp_count + del_count)) >= filter_thresh:
            remove_opposite_cnas(filtered_cnas, gene_cnas, DEL)
        else:
            for cna in gene_cnas:
                filtered_cnas.remove(cna)
                
    return filtered_cnas

def remove_opposite_cnas(cnas, gene_cnas, mut_type):
    for cna in gene_cnas:
        if cna.mut_type != mut_type:
            cnas.remove(cna)

def mut_heat(num_samples, snvs, cnas, min_freq):
    """Return a dict mapping gene name to heat score based on the given mutation data.
    
    Arguments:
    num_samples -- the number of samples tested for mutations
    snvs -- iterable of Mutation tuples representing SNVs (mut_type == SNV)
    cnas -- iterable of Mutation tuples representing CNAs (mut_type == AMP or DEL)
    min_freq -- the minimum number of samples in which a gene must have an SNV to be considered
                mutated in the heat score calculation.
    
    """
    
    genes2mutations = defaultdict(set)
    for snv in snvs:
        genes2mutations[snv.gene].add(snv)
    for cna in cnas:
        genes2mutations[cna.gene].add(cna)
    
    print("\t- Including %s genes in %s samples at min frequency %s" %
          (len(genes2mutations), num_samples, min_freq))
    
    return dict([(g, len(heat) / float(num_samples)) for g, heat in genes2mutations.items()
                 if num_snvs(heat) >= min_freq or num_cnas(heat) > 0])

def expr_filter_heat(gene2heat, genes_to_preserve):
    """Return (1) a dict mapping gene names to heat scores containing only entries for genes in
    genes_to_preserve, and (2) a list of genes whose heat scores are not included in the returned dict.
    
    Arguments:
    gene2heat -- dict mapping gene names to heat scores
    genes_to_preserve -- set of genes whose heat scores should be contained in the returned dict
    
    """
    filtered_heat = dict()
    excluded_genes = list()
    for gene in gene2heat:
        if gene in genes_to_preserve:
            filtered_heat[gene] = gene2heat[gene]
        else:
            excluded_genes.append(gene)
    return filtered_heat, excluded_genes
