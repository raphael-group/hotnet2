from collections import defaultdict, namedtuple
import multiprocessing as mp
import random
from constants import Mutation, SNV
import heat
import hnio

################################################################################
# Heat permutation

def heat_permutation_wrapper((heat_scores, eligible_genes)):
    permuted_genes = list(eligible_genes)
    random.shuffle(permuted_genes)
    permuted_genes = permuted_genes[:len(heat_scores)]

    permuted_heat = dict((gene, heat) for gene, heat in zip(permuted_genes, heat_scores))

    return permuted_heat

def permute_heat(heat, network_genes, num_permutations, addtl_genes=None, num_cores=1):
    """Return a list of num_permutation dicts, each mapping gene names to permuted heat scores.
    
    Arguments:
    heat -- dict mapping a gene name to the heat score for that gene
    network_genes -- iterable of names of genes in the network
    num_permutations -- number of heat permutations to produce
    addtl_genes -- iterable of names of genes that do not have heat scores in the real data but
                   which may have heat scores assigned in permutations. Defaults to None.
    num_cores -- number of cores to use for running in parallel
    
    """
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    heat_scores = heat.values()
    if not addtl_genes: addtl_genes = set()
    genes_eligible_for_heat = set(heat.keys()).union(addtl_genes).intersection(network_genes)
    
    args = [(heat_scores, genes_eligible_for_heat)] * num_permutations
    permutations = map_fn(heat_permutation_wrapper, args)
    
    if num_cores != 1:
        pool.close()
        pool.join()

    return permutations

################################################################################
# Mutation permutation

def mutation_permuation_heat_wrapper((samples, genes, cnas, gene2length, bmr, gene2bmr, gene2chromo,
                                     chromo2genes, cna_filter_threshold, min_freq)):
    permuted_snvs = permute_snvs(samples, genes, gene2length, bmr, gene2bmr)
    permuted_cnas = permute_cnas(cnas, gene2chromo, chromo2genes)
    if cna_filter_threshold:
        permuted_cnas = heat.filter_cnas(permuted_cnas, cna_filter_threshold)
        
    return heat.mut_heat(genes, len(samples), permuted_snvs, permuted_cnas, min_freq)
    

def generate_mutation_permutation_heat(heat_fn, sample_file, gene_file, genes_in_network, snv_file,
                                       gene_length_file, bmr, bmr_file, cna_file, gene_order_file,
                                       cna_filter_threshold, min_freq, num_permutations, num_cores=1):
    """Return a list of num_permutation dicts, each mapping gene names to heat scores calculated
    from permuted mutation data.
    
    Arguments:
    heat_fn -- the name of the function used to calculate heat scores for the real data.
               A RuntimeError is raised if this is not "load_mutation_heat"
    sample_file -- path to TSV file containing tested sample IDs as the first column. If None, the
                   set of samples is assumed to be all samples that are provided in the SNV or CNA data.
    gene_file -- path to file containing names of tested genes, one per line. If None, the set of
                 tested genes is assumed to be all genes that have mutations in either the SNV or
                 CNA data.
    genes_in_network -- iterable of names of genes in the PPI network
    snv_file -- path to TSV file containing SNVs where the first column of each line is a sample ID
                and subsequent columns contain the names of SNVs with mutations in that sample.
    gene_length_file -- path to TSV file containing gene names in the first column and the length
                        of the gene in base pairs in the second column
    bmr -- default background mutation rate
    bmr_file -- path to TSV file with gene names in the first column and the background mutation rate
                for the gene in the second column. The default BMR will be used for any gene without
                a BMR in this file. If None, the default BMR will be used for all genes.
    cna_file -- path to TSV file containing CNAs where the first column of each line is a sample ID
                and subsequent columns contain gene names followed by "(A)" or "(D)" indicating an
                ammplification or deletion in that gene for the sample. Lines starting with '#'
                will be ignored.
    gene_order_file -- path to file containing tab-separated lists of genes on each chromosme,
                       one chromosome per line
    cna_filter_threshold -- proportion of CNAs in a gene across samples that must share  the same
                            CNA type in order for the CNAs to be included. Must be > .5 if not None.
                            If None, all CNAs will be included.
    min_freq -- the minimum number of samples in which a gene must have an SNV to be considered
                mutated in the heat score calculation.
    num_permutations -- the number of permuted data sets to be generated
    num_cores -- number of cores to use for running in parallel
    
    """
    if heat_fn != "load_mutation_heat":
        raise RuntimeError("Heat scores must be based on mutation data to perform\
                            delta selection based on mutation data permutation.")
    
    samples = hnio.load_samples(sample_file) if sample_file else None
    genes = hnio.load_genes(gene_file) if gene_file else None
    cnas = hnio.load_cnas(cna_file, genes, samples)
    gene2length = hnio.load_gene_lengths(gene_length_file)
    gene2bmr = hnio.load_gene_specific_bmrs(bmr_file) if bmr_file else {}
    gene2chromo, chromo2genes = hnio.load_gene_order(gene_order_file)
    
    if not samples:
        snvs = hnio.load_snvs(snv_file, genes, samples)
        samples = set([snv.sample for snv in snvs] + [cna.sample for cna in cnas])
    if not genes:
        genes = set([snv.gene for snv in snvs] + [cna.gene for cna in cnas])
    
    #only generate mutations for genes that are in the network
    genes = set(genes).intersection(genes_in_network)
    
    if num_cores != 1:
        pool = mp.Pool(None if num_cores == -1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map
    
    args = [(samples, genes, cnas, gene2length, bmr, gene2bmr, gene2chromo,chromo2genes,
             cna_filter_threshold, min_freq)] * num_permutations
    heat_permutations = map_fn(mutation_permuation_heat_wrapper, args)
    
    if num_cores != 1:
        pool.close()
        pool.join()
    
    return heat_permutations

def permute_snvs(samples, tested_genes, gene2length, bmr, gene2bmr):
    permuted_snvs = []
    for sample in samples:
        for gene in tested_genes:
            gene_bmr = gene2bmr[gene] if gene in gene2bmr else bmr
            gene_length = gene2length[gene]
            prob = 1 - pow(1 - gene_bmr, gene_length)
            if random.random() <= prob:
                permuted_snvs.append(Mutation(sample, gene, SNV))

    return permuted_snvs

Block = namedtuple("Block", ["chromosome", "start_index", "mut_type", "genes"])

def permute_cnas(cnas, gene2chromo, chromo2genes):
    samples2cnas = defaultdict(list)
    for cna in cnas:
        samples2cnas[cna.sample].append(cna)
    
    permuted_cnas = []
    for sample in samples2cnas:
        chromo2blocks = get_cna_blocks_for_sample(samples2cnas[sample], gene2chromo, chromo2genes)
        for chromo, blocks in chromo2blocks.iteritems():
            genes = chromo2genes[chromo]
            invalid_indices = []
            for block in blocks:
                permuted_indices = get_block_indices(len(genes), len(block.genes), invalid_indices)
                for index in permuted_indices:
                    permuted_cnas.append(Mutation(sample, genes[index], block.mut_type))

                new_invalid_indices = permuted_indices +\
                    [min(permuted_indices) - 1, max(permuted_indices) + 1]
                invalid_indices.extend(new_invalid_indices)

    return permuted_cnas

def get_block_indices(chromo_length, block_length, invalid_indices, max_attempts=10):
    start = random.randint(0, chromo_length - block_length)
    attempts = 1
    while not is_block_valid(start, block_length, invalid_indices):
        if attempts == max_attempts:
            raise RuntimeError("Cannot place CNA block after %s attempts" % (max_attempts))
        start = random.randint(0, chromo_length - block_length)
        attempts += 1

    return range(start, start + block_length)
    

def is_block_valid(start, block_length, invalid_indices):
    for i in range(start, start + block_length):
        if i in invalid_indices:
            return False
    return True

def get_cna_blocks_for_sample(cnas, gene2chromo, chromo2genes):
    #sort cnas by chromosome, then by position in chromosome
    cnas.sort(key = lambda cna: (gene2chromo[cna.gene],
                                 chromo2genes[gene2chromo[cna.gene]].index(cna.gene)))
    
    chromo2blocks = defaultdict(list)
    curr_chromo = gene2chromo[cnas[0].gene]
    curr_start_index = chromo2genes[curr_chromo].index(cnas[0].gene)
    curr_mut_type = cnas[0].mut_type
    curr_genes = [cnas[0].gene]
    for cna in cnas[1:]:
        cna_index = chromo2genes[gene2chromo[cna.gene]].index(cna.gene)
        if (gene2chromo[cna.gene] == curr_chromo and
            cna_index == curr_start_index + len(curr_genes) and
            cna.mut_type == curr_mut_type):
            curr_genes.append(cna.gene)
        else:
            #store completed block
            chromo2blocks[curr_chromo].append(Block(curr_chromo, curr_start_index,
                                                    curr_mut_type, curr_genes))

            #start new block
            curr_chromo = gene2chromo[cna.gene]
            curr_start_index = cna_index
            curr_mut_type = cna.mut_type
            curr_genes = [cna.gene]
    chromo2blocks[curr_chromo].append(Block(curr_chromo, curr_start_index,
                                            curr_mut_type, curr_genes))

    return chromo2blocks
