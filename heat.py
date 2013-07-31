from math import log10
from collections import defaultdict
import scipy
from constants import *

def num_snvs(mutations):
    return len([mut for mut in mutations if mut.mut_type == SNV])

def num_cnas(mutations):
    return len([mut for mut in mutations if mut.mut_type == AMP or mut.mut_type == DEL])

def filter_cnas(cnas, filter_thresh):
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
    genes2mutations = defaultdict(set)
    for snv in snvs:
        genes2mutations[snv.gene].add(snv)
    for cna in cnas:
        genes2mutations[cna.gene].add(cna)
    
    print("\t- Including %s genes in %s samples at min frequency %s" %
          (len(genes2mutations), num_samples, min_freq))
    
    return dict([(g, len(heat) / float(num_samples)) for g, heat in genes2mutations.items()
                 if num_snvs(heat) >= min_freq or num_cnas(heat) > 0])

def mut_heat_old(samples2snvs, samples2cnas, min_freq):
    n = float(len(set(samples2snvs.keys()) | set(samples2cnas.keys())))
    
    genes2mutations = defaultdict(set)
    snvs = reduce(lambda acc, upd: acc.union(upd), samples2snvs.values(), set())
    cnas = reduce(lambda acc, upd: acc.union(upd), samples2cnas.values(), set())
    for snv in snvs:
        genes2mutations[snv.gene].add(snv)
    for cna in cnas:
        genes2mutations[cna.gene].add(cna)
    
    print "\t- Including %s genes in %s samples at min frequency %s" % (len(genes2mutations), int(n), min_freq)
    
    return dict([(g, len( heat ) / n) for g, heat in genes2mutations.items()
                 if num_snvs(heat) >= min_freq or num_cnas(heat) > 0])

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
    gene2heat = dict([(gene, -log10(score["qval"]))
                      for gene, score in gene2mutsig.items()
                      if score["qval"] < threshold])
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat

def music_heat(gene2music, threshold=1.0, max_heat=15):
    print "* Creating MuSiC heat map..."
    print "\t- Mapping q-values of 0 to", max_heat
    def music_heat(qvals):
        heat = scipy.median([ qvals["FDR_CT"], qvals["FDR_LRT"], qvals["FDR_FCPT"] ])
        return -log10(heat) if heat != 0 else max_heat
    gene2heat = dict([(gene, music_heat(scores)) for gene, scores in gene2music.items()
                      if scipy.median(scores.values()) < threshold])
    print "\t- Including", len(gene2heat), "genes at threshold", threshold
    return gene2heat

def expr_filter_heat(gene2heat, genes_to_preserve):
    gene2heat = dict([(g, h) for g, h in gene2heat.items() if g in genes_to_preserve])
    return gene2heat
