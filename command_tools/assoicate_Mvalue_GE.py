#!/usr/bin/env python3
# coding=utf-8
"""Correlate peak M-values with target-gene expression fold changes."""
import argparse

from scipy.stats import pearsonr

from import_MAmotif import GeneSet, read_refgenes, read_MAnorm_peaks
from io_related import read_fold_change_file


def cal_pearson_correlation(pk_fp, refgene_fp, genes_fold_change_fp):
    """Pearson correlation between each peak M-value and its target gene log2FC."""
    ref_genes = GeneSet(read_refgenes(refgene_fp))
    pks = read_MAnorm_peaks(pk_fp)
    targenes = [ref_genes.find_target_gene(pk) for pk in pks]
    pks = [pks[i] for i in range(len(targenes)) if len(targenes[i]) > 0]
    targenes = [genes[0] for genes in targenes if len(genes) > 0]

    gene_log2fcs = read_fold_change_file(genes_fold_change_fp)
    m_values = [pk.mvalue for pk in pks]
    log2fcs = [
        gene_log2fcs[gene.name2] if gene.name2 in gene_log2fcs else None
        for gene in targenes
    ]
    m_values = [m_values[i] for i in range(len(pks)) if log2fcs[i] is not None]
    log2fcs = [value for value in log2fcs if value is not None]

    prs = pearsonr(m_values, log2fcs)
    print(prs)

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(m_values, log2fcs)
    plt.xlabel('mvalue')
    plt.ylabel('log2 fold change')
    plt.title('pearsonr: %f' % prs[0])
    plt.savefig(pk_fp[:-4] + '.png')


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='MAnorm peak file path')
    parser.add_argument('-r', dest='refgene', required=True, help='refgene file')
    parser.add_argument('-f', dest='fold_change', required=True, help='log2 fold-change file')
    options = parser.parse_args(argv)
    cal_pearson_correlation(options.pk, options.refgene, options.fold_change)
    print('Done!')


if __name__ == '__main__':
    command()
