from collections import defaultdict
from math import log10
import scipy
from constants import Mutation, SNV, AMP, DEL

def filter_heat(heat, min_score, zero_genes=False, message=None):
    """Returns (1) a dict mapping gene names to heat scores where all non-zero heat scores are at
    least min_score and any genes that originally had score less than min_score either have score 0
    or are excluded depending on the value of zero_genes; and (2) a set of the genes that had scores
    less than min_score

    Arguments:
    heat -- dict mapping gene names to heat scores
    min_score -- minimum heat score a gene must have to be included or non-zero (depending on the
                 value of zero_genes) in the returned dict. If None, the minimum score will be
                 calculated as the minimum of the non-zero heat scores included in the input dict
    zero_genes -- if true, genes with score below min_score will be included in the returned heat
                  dict with their scores set to zero; otherwise, they will be excluded
    message -- message to print about number of filtered genes; '##' will be replaced with the number

    """
    if not min_score:
        min_score = min([score for score in heat.values() if score > 0])

    filtered_genes = set()
    filtered_heat = dict()
    for gene, score in heat.iteritems():
        if score >= min_score:
            filtered_heat[gene] = score
        else:
            filtered_genes.add(gene)
            if zero_genes:
                filtered_heat[gene] = 0

    if message and len(filtered_genes) > 0:
        print '\t- ' + message.replace('##', str(len(filtered_genes)))

    return filtered_heat, filtered_genes

def num_snvs(mutations):
    """Return the number of valid SNVs in the given iterable of Mutations
    (those with mut_type == SNV and valid == True).

    Arguments:
    mutations -- iterable of Mutation tuples

    """
    return len([mut for mut in mutations if mut.mut_type == SNV and mut.valid])

def num_cnas(mutations):
    """Return the number of valid CNAs in the given iterable of Mutations
    (those with mut_type == AMP or DEL and valid == True).

    Arguments:
    mutations -- iterable of Mutation tuples

    """
    return len([mut for mut in mutations if (mut.mut_type == AMP or mut.mut_type == DEL) and mut.valid])

def filter_cnas(cnas, filter_thresh):
    """Return a list of Mutation tuples in which CNAs in genes where the proportion of CNAs in the
    gene across samples is less than filter_thresh have their "valid" field set to False. CNAs in
    genes whose CNAs pass the threshold that are of the opposite type of the dominant type in the gene
    also have their "valid" field set to False.

    Arguments:
    cnas -- a list of Mutation tuples representing CNAs. All are assumed to be of mut_type AMP or DEL
    filter_thresh -- the proportion of CNAs in a gene that must be of the same type

    """
    if filter_thresh <= .5:
        raise ValueError("filter_thresh must be greater than .5")

    filtered_cnas = list()
    genes2cnas = defaultdict(list)
    for cna in cnas:
        genes2cnas[cna.gene].append(cna)

    for gene_cnas in genes2cnas.itervalues():
        amp_count = float(len([cna for cna in gene_cnas if cna.mut_type == AMP]))
        del_count = float(len([cna for cna in gene_cnas if cna.mut_type == DEL]))
        if (amp_count / (amp_count + del_count)) >= filter_thresh:
            invalidate_opposite_cnas(filtered_cnas, gene_cnas, AMP)
        elif (del_count / (amp_count + del_count)) >= filter_thresh:
            invalidate_opposite_cnas(filtered_cnas, gene_cnas, DEL)
        else:
            for cna in gene_cnas:
                filtered_cnas.append(get_invalidated_mutation(cna))

    return filtered_cnas

def invalidate_opposite_cnas(cnas, gene_cnas, mut_type):
    for cna in gene_cnas:
        if cna.mut_type == mut_type:
            cnas.append(cna)
        else:
            cnas.append(get_invalidated_mutation(cna))

def get_invalidated_mutation(mutation):
    return Mutation(mutation.sample, mutation.gene, mutation.mut_type, False)

def mut_heat(genes, num_samples, snvs, cnas, min_freq):
    """Return a dict mapping gene name to heat score based on the given mutation data.

    Arguments:
    genes -- iterable of genes tested for mutations
    num_samples -- the number of samples tested for mutations
    snvs -- iterable of Mutation tuples representing SNVs (mut_type == SNV)
    cnas -- iterable of Mutation tuples representing CNAs (mut_type == AMP or DEL)
    min_freq -- the minimum number of samples in which a gene must have an SNV to be considered
                mutated in the heat score calculation.

    """

    genes2mutations = dict((gene, set()) for gene in genes)
    for snv in snvs:
        genes2mutations[snv.gene].add(snv)
    for cna in cnas:
        genes2mutations[cna.gene].add(cna)
    print("* Calculating heat scores for %s genes in %s samples at min frequency %s" %
          (len(genes2mutations), num_samples, min_freq))

    gene2heat = dict()
    for gene, mutations in genes2mutations.iteritems():
        snv_mut_samples = set( m.sample for m in mutations if m.mut_type == SNV and m.valid )
        cna_mut_samples = set( m.sample for m in mutations if (m.mut_type == AMP or m.mut_type == DEL) and m.valid )

        # Minimum frequency is for SNVs *only*, so we just use CNAs if the SNVs
        # are below min_freq
        if len(snv_mut_samples) < min_freq:
            gene2heat[gene] = len(cna_mut_samples) / float(num_samples)
        else:
            gene2heat[gene] = len(snv_mut_samples | cna_mut_samples) / float(num_samples)

    return gene2heat

NULL = 100
def fm_heat(gene2heat, fm_threshold, cis_threshold=0.01, CIS=False):
    print "* Creating oncodrive heat map..."
    if CIS: print "\tIncluding CIS scores at threshold", cis_threshold, "..."
    heat = dict()
    src_fm, src_cis_amp, src_cis_del = 0, 0, 0
    for g, scores in gene2heat.items():
        if CIS:
            del_score = scores["del"] if scores["del"] < cis_threshold else NULL
            amp_score = scores["amp"] if scores["amp"] < cis_threshold else NULL
            fm_score  = scores["fm"] if scores["fm"] < fm_threshold else NULL
            if fm_score == NULL and amp_score == NULL and del_score == NULL: continue
            min_val = min(del_score, amp_score, fm_score)
            heat[g] = -log10( min_val )
            if min_val == scores["fm"]: src_fm += 1
            elif min_val == scores["amp"]: src_cis_amp += 1
            elif min_val == scores["del"]: src_cis_del += 1
        else:
            if scores["fm"] >= fm_threshold: continue
            heat[g] = -log10(scores["fm"])
            src_fm += 1
    print "\t- Genes using FM score:", src_fm
    print "\t- Genes using CIS AMP score:", src_cis_amp
    print "\t- Genes using CIS DEL score:", src_cis_del

    return heat

def mutsig_heat(gene2mutsig, threshold=1.0):
    print "* Creating MutSig heat map..."
    gene2heat = dict((gene, -log10(score["qval"]))
                     for gene, score in gene2mutsig.items()
                     if score["qval"] < threshold)
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat

def music_heat(gene2music, threshold=1.0, max_heat=15):
    print "* Creating MuSiC heat map..."
    print "\t- Mapping q-values of 0 to", max_heat
    def music_heat(qvals):
        heat = scipy.median([ qvals["FDR_CT"], qvals["FDR_LRT"], qvals["FDR_FCPT"] ])
        return -log10(heat) if heat != 0 else max_heat
    gene2heat = dict((gene, music_heat(scores)) for gene, scores in gene2music.items()
                     if scipy.median(scores.values()) < threshold)
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat

def filter_heat_to_network_genes(gene2heat, network_genes, verbose):
    """Return a dict mapping genes to heat scores such that only genes in the network are included.

    Arguments:
    gene2heat -- dict mapping gene names to heat scores
    network_genes -- set of network_genes

    """
    filtered_heat = dict()
    num_removed = 0
    for gene, heat in gene2heat.iteritems():
        if gene in network_genes:
            filtered_heat[gene] = heat
        else:
            num_removed += 1

    if verbose > 1 and num_removed > 0:
        print "\t- Removing %s genes not in the network" % num_removed

    return filtered_heat

def reconcile_heat_with_tested_genes(gene2heat, tested_genes):
    """Return a dict mapping gene names to heat scores containing for each gene in tested_genes
    and only for genes in tested_genes. Genes in tested_genes not in gene2heat will be given a
    score of 0.

    Arguments:
    gene2heat -- dict mapping gene names to heat scores
    tested_genes -- set of genes that should have heat scores in the returned dict

    """
    filtered_heat = dict((g, gene2heat[g] if g in gene2heat else 0) for g in tested_genes)

    num_removed = len(set(gene2heat.keys()).difference(tested_genes))
    if num_removed > 0:
        print "\t- Removing %s genes not in gene_filter_file" % num_removed

    num_zeroed = len(set(tested_genes).difference(gene2heat.keys()))
    if num_zeroed > 0:
        print "\t- Assigned score 0 to %s genes in gene file without scores" % num_zeroed

    return filtered_heat
