# coding=utf-8
import os
from math import log10
from import_MAmotif import read_MAnorm_peaks, read_refgenes, \
    GeneSet, make_2geneset_fisher_exact_test, read_gene_list_file
import matplotlib.pyplot as plt
import numpy as np


__author__ = 'Semal'


def get_specific_line(fp, string='ZNF263'):
    with open(fp) as fi:
        for line in fi:
            if string in line:
                return line

# 抽取Hek293和其他细胞系比较的MAmotif结果
def extract_res():
    os.chdir('F:\\2.MAmotif\\10. HEK293\\4.MAmotif\H3K4me3_Hek293_VS_otherCellLines')
    for fd_name in os.listdir('.'):
        for fn in os.listdir(fd_name):
            if 'Hek293' not in fn or 'MAmotif' not in fn:
                continue
            fp = os.sep.join([fd_name, fn])
            line = get_specific_line(fp)
            print fp, ':', line,


# 将抽取的结果按照distal,promoter,both分开
def divide_3part():
    os.chdir('F:\\2.MAmotif\\10. HEK293\\4.MAmotif')
    distal_res, promoter_res, both_res = '', '', ''
    with open('H3K4me3_Hek293_ZNF263_MAmotif.txt') as fi:
        for line in fi:
            if 'distal' in line:
                distal_res += line
                continue
            if 'promoter' in line:
                promoter_res += line
                continue
            both_res += line
    open('H3K4me3_Hek293_ZNF263_MAmotif_distal.txt', 'w').write(distal_res)
    open('H3K4me3_Hek293_ZNF263_MAmotif_promoter.txt', 'w').write(promoter_res)
    open('H3K4me3_Hek293_ZNF263_MAmotif_both.txt', 'w').write(both_res)


# 将上面的结果处理成能够用来做图表的格式
def process_res():
    os.chdir('F:\\2.MAmotif\\10. HEK293\\4.MAmotif')
    with open('H3K4me3_Hek293_ZNF263_MAmotif_promoter.txt') as fi:
        for line in fi:
            cnt = line.split('\t')
            print cnt[0].replace('H3K4me3_Hek293_UW_', ''), '\t',cnt[-1],


# 作图表
def show_res():
    his_mod = 'H3K9ac'
    motif = 'MZF1_5-13'
    os.chdir('F:\\2.MAmotif\\2. MAnorm2\\H1hESC_%s' % his_mod)
    # print 'Current directory:', 'F:\\2.MAmotif\\2. MAnorm2\\Esb4_%s' % his_mod
    import re
    pattern = r'[a-zA-Z0-9]*_vs_[a-zA-Z0-9]*'
    def get_info(fn, promoter_distal='promoter'):
        names, values = [], []
        with open(fn) as fi:
            for line in fi:
                if line.startswith('#'):
                    continue
                if promoter_distal in line and motif in line:
                    cnt = line.split('\t')
                    name = re.search(pattern, cnt[0]).group()
                    names.append(name.replace('_', ' ').replace('promoter', '').replace('distal', ''))
                    values.append(-log10(float(cnt[-1].strip())))
        return names, values

    names, values = get_info('4motif_MAmotif.txt')
    new_zipper = sorted(zip(values, names))
    names = [z[1] for z in new_zipper]
    values = [z[0] for z in new_zipper]

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(131)
    ax.barh(range(1, 2 * len(values) + 1, 2), values, 1.5, align='center', color='b')
    ax.set_yticks(range(1, 2 * len(values) + 1, 2))
    ax.set_yticklabels(names, fontsize=15)
    ax.set_xlim(0, 10)
    ax.set_xticks([0, 5, 10])
    ax.set_title('Promoter %s' % his_mod, fontdict={'fontsize': 18})
    # ax.tight_layout()
    ax.set_xlabel('-log10(P-Value of %s)' % motif, fontsize=15)
    # ax.ylabel('Comparisons')
    ax.plot([3, 3], [0, 2 * len(values)], '--', color='r')
    ax.text(3, 0.5, '3')
    # ax.show()
    # fig_name = raw_input('Name this figure!\n')
    # plt.savefig(fig_name)

    names, values = get_info('4motif_MAmotif.txt', 'distal')
    new_zipper = sorted(zip(values, names))
    names = [z[1] for z in new_zipper]
    values = [z[0] for z in new_zipper]

    ax = fig.add_subplot(133)
    ax.barh(range(1, 2 * len(values) + 1, 2), values, 1.5, align='center', color='b')
    ax.set_yticks(range(1, 2 * len(values) + 1, 2))
    ax.set_yticklabels(names, fontsize=15)
    ax.set_xlim(0, 10)
    ax.set_xticks([0, 5, 10])
    ax.set_title('Distal %s' % his_mod, fontdict={'fontsize': 18})
    # ax.tight_layout()
    ax.set_xlabel('-log10(P-Value of %s)' % motif, fontsize=15)
    # ax.ylabel('Comparisons')
    ax.plot([3, 3], [0, 2 * len(values)], '--', color='r')
    ax.text(3, 0.5, '3')
    plt.show()
    # fig_name = raw_input('Name this figure!\n')
    # ax.savefig(fig_name)
    # plt.savefig('H1hESC_ZNF263.png')


# -------------------------------- GSEA-like -------------------------------
def connect_M_and_genes(**kwargs):
    peak = kwargs.pop('peak')
    refgene = kwargs.pop('refgene')

    peaks = read_MAnorm_peaks(peak)
    refgenes = GeneSet(read_refgenes(refgene))

    M_gene_pairs = []
    for peak in peaks:
        target = refgenes.find_target_gene(peak, True)
        if len(target) == 0:
            continue
        M_gene_pairs.append((peak.mvalue, target[0].name2))

    return M_gene_pairs


def scan_and_calculate(M_gene_pairs, DE_genes, **kwargs):

    M_gene_pairs = sorted(M_gene_pairs)[::-1]
    M = [pair[0] for pair in M_gene_pairs]
    genes = [pair[1] for pair in M_gene_pairs]
    background_genes = set(genes)

    window = kwargs.pop('window', 100)
    step = kwargs.pop('step', 1)

    enrich_score = []
    for i in xrange(0, len(M) - window, step):
        window_genes = genes[i: i + window]
        es = make_2geneset_fisher_exact_test(background_genes, set(window_genes), set(DE_genes))
        enrich_score.append(es)

    # plot enrich score curve
    plt.figure(21)
    ax1 = plt.subplot(211)
    ax1.plot(range(len(enrich_score)), enrich_score, color='g')
    ax1.plot([0, len(enrich_score)], [1, 1], '--')
    ax1.set_ylabel('enrichment score')
    # plot M value heat bar
    # ax2 = plt.subplot(312)
    # a = range(0, len(enrich_score))
    # b = [0, 1]
    # from pylab import meshgrid
    # x, y = meshgrid(a, b)
    # ax2.pcolor(x, y, np.array([M[:len(enrich_score)], M[:len(enrich_score)]]))
    # plot M value curve
    ax2 = plt.subplot(212)
    ax2.plot(range(0, len(enrich_score)), M[:len(enrich_score)], color='grey')
    ax2.plot([0, len(enrich_score)], [0, 0], '--')
    ax2.set_ylabel('M')
    ax2.set_xlabel('M-ranked-peak')
    plt.show()


def test():
    peak_file_path = 'F:\\2.MAmotif\\10. ' \
                     'HEK293\\2.MAnorm\\H3K4me3_Hek293_UW_vs_H1hESC_Braod_Rep1_P100' \
                     '\\Hek293_H3K4me3_UW_Rep1_P100_peaks_MAvalues.xls'
    refgene_file_path = 'F:\\MAmotif_src\\data\\refGene_hg19.txt'
    DE_gene = \
        read_gene_list_file('F:\\2.MAmotif\\10. '
                            'HEK293\\5.Hek293_RNAseq\\HEK293_mRNA_UMS_H1hesc_DEseq_DE_log2FC_-1_p0'
                            '.01_geneSymbols.txt')
    paris = connect_M_and_genes(peak=peak_file_path, refgene=refgene_file_path)
    scan_and_calculate(paris, DE_gene, window=1000)


if __name__ == '__main__':
    # extract_res()
    # divide_3part()
    show_res()