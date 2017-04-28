# coding=utf-8
# 我们将peaks按照是否有指定motif分成两组，然后将两组peaks的靶基因找出来，对这两组的靶基因的log2(fold_change)值做
# 分布上的差异检验（类似于两组M值间的分布差异检验）
import os
from scipy import stats
import numpy as np

from import_MAmotif import GeneSet, read_refgenes, match_manorm_with_motifscan
from io_related import read_fold_change_file


def do_ttest(yes, no):
        try:
            t_statistic, two_tailed_pvalue = stats.ttest_ind(yes, no, equal_var=False)
            if t_statistic < 0:
                pvalue_left = two_tailed_pvalue / 2
                pvalue_right = 1 - two_tailed_pvalue / 2
            else:
                pvalue_left = 1 - two_tailed_pvalue / 2
                pvalue_right = two_tailed_pvalue / 2
            return t_statistic, pvalue_left, pvalue_right
        except:
            return None


def do_ranksum_test(yes, no):
    try:
        z_statistic, two_pvalue = stats.ranksums(yes, no)
        if z_statistic < 0:
            pvalue_left = two_pvalue / 2
            pvalue_right = 1 - two_pvalue / 2
        else:
            pvalue_left = 1 - two_pvalue / 2
            pvalue_right = two_pvalue / 2
        return z_statistic, pvalue_left, pvalue_right
    except:
        return None


def motifs_classified_peaks_targenes_ttest(peak_fp, refgene_fp, genes_fc_fp, motifscan_fp):
    """
    this method will match peaks with its target genes, then using motif_scan_result to classify peaks. the test result
    of corresponded target genes will be output in a table.
    @param peak_fp: a peak file
    @param genes_fc_fp: the 2 cells of this fold change file should be the same of peak file
    """
    pks, tarnum_dict, motifs = match_manorm_with_motifscan(peak_fp, motifscan_fp)
    refgenes = GeneSet(read_refgenes(refgene_fp))
    targenes = [refgenes.find_target_gene(pk) for pk in pks]

    log2fc = read_fold_change_file(genes_fc_fp)

    # classify pks
    pk_name = os.path.split(peak_fp)[1]
    test_result = open(pk_name[:-4] + '_motif_classified_target_genes_test_result' + '.xls', 'w')
    header = \
        'Motif Name\t' \
        'Target Number\tAverage of Target log2-fold-change\tDeviation of Target log2-fold-change\t' \
        'Non-target Number\tAverage of Non-target log2-fold-change\tDeviation of Non-target log2-fold-change\t' \
        'T-test Statistics\tT-test P-value(right-tail)\t' \
        'RanSum-test Statistics\tRankSum-test P-value(right-tail)\t' \
        'Maximal P-value(T-test,RankSum-test)\t' \
        'P-value By Benjamin correction\tP-value By Bonferroni-correction\n'
    test_result.write(header)
    lines = []
    ttest_result, rtest_result = [], []
    max_pvalue = []

    def get_genes(targenes_index_array):
        genes = []
        for index in targenes_index_array:
            for gene in targenes[index]:
                genes.append(gene.name2)
        return set(genes)

    def get_genes_log2fc(genes_name):
        return np.array([log2fc[name] for name in genes_name if name in log2fc.keys()])

    for moti in motifs:
        yes_gene_names = get_genes(np.where(tarnum_dict[moti] > 0)[0])
        yes = get_genes_log2fc(yes_gene_names)
        no_gene_names = get_genes(np.where(tarnum_dict[moti] == 0)[0])
        no = get_genes_log2fc(no_gene_names)

        line = '%s\t' % moti
        line += '{0:d}\t{1:f}\t{2:f}\t'.format(yes.size, yes.mean(), yes.std())
        line += '{0:d}\t{1:f}\t{2:f}\t'.format(no.size, no.mean(), no.std())
        lines.append(line)
        print line
        ttest = do_ttest(yes, no)
        ttest_result.append(ttest)
        rtest = do_ranksum_test(yes, no)
        rtest_result.append(rtest)
        if ttest is not None:
            max_pvalue.append(max(ttest[2], rtest[2]))
        else:
            max_pvalue.append(None)

    idx_sorted = np.array(max_pvalue).argsort().tolist()
    for i in idx_sorted:
        line = lines[i]
        if ttest_result[i] is not None:
            line += str(ttest_result[i][0]) + '\t' + str(ttest_result[i][2]) + '\t'
            line += str(rtest_result[i][0]) + '\t' + str(rtest_result[i][2]) + '\t'
            line += \
                str(max_pvalue[i]) + '\t' + \
                str(len(motifs) * max_pvalue[i] / (i + 1) if len(motifs) * max_pvalue[i] / (i + 1) < 1 else 1) + '\t' +\
                str(len(motifs) * max_pvalue[i] if len(motifs) * max_pvalue[i] < 1 else 1) + '\n'
        else:
            line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (None, None, None, None, None, None, None)
        # print line
        test_result.write(line)


def command():
    from optparse import OptionParser

    opt_parser = OptionParser()
    opt_parser.add_option('-p', dest='pk', help='cell to cell peak file path')
    opt_parser.add_option('-r', dest='refgene', help='refgene file')
    opt_parser.add_option('-f', dest='fold_change', help='cell to cell log2 fold change')
    opt_parser.add_option('-m', dest='motif', help='motifscan result file path')
    options, args = opt_parser.parse_args()
    pk = options.pk
    refgene = options.refgene
    fold_change = options.fold_change
    motif = options.motif
    motifs_classified_peaks_targenes_ttest(pk, refgene, fold_change, motif)
    print 'Done!'


def test_motifs_classified_peaks_targenes_ttest():
    from constant import manorm_peak_file as pk
    from constant import hg19_refgenes_file as ref
    from constant import gene_fold_change_file as fc
    from constant import motifscan_file as ms
    motifs_classified_peaks_targenes_ttest(pk, ref, fc, ms)


if __name__ == '__main__':
    # command()
    test_motifs_classified_peaks_targenes_ttest()
    pass