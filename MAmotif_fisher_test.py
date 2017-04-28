# coding=utf-8
# 富集分析相关
from scipy.stats.stats import fisher_exact
import os
from MAmotif_pkg import *
from gene_set_analysis import make_2geneset_fisher_exact_test, get_fisher_table, cal_enrich_score
from io_related import read_gene_list_file


def get_motifs_target_genes(bkg_peak_file, motifscan_result, refgene_file, is_nearest=False):
    """
    获取作为背景的peak的靶基因symbol，以及motifscan中所有motif的靶基因symbol
    @param bkg_peak_file: 用来产生背景基因的MAnorm peak文件
    @param motifscan_result: motifscan的tarnum文件
    @param refgene_file: refseq基因文件
    @return: 背景基因和所有motif的靶基因symbol组成的字典
    """
    # get pk list, target gene list and target number list
    pk_list, target_num_dict, motifs = match_manorm_with_motifscan(bkg_peak_file, motifscan_result)
    gene_list = read_refgenes(refgene_file)
    gene_set = GeneSet()
    gene_set.set_sequences(gene_list)
    print 'get peaks target genes...'
    target_genes_list = [gene_set.find_target_gene(pk, is_nearest) for pk in pk_list]

    all_bk_genes = []
    for ge in target_genes_list:
        all_bk_genes += ge
    bk_genes = [gene.name2 for gene in all_bk_genes]
    bk_genes = set(bk_genes)

    motifs_genes = {}
    for motif in motifs:
        other_genes = []
        for tarnum, ge in zip(target_num_dict[motif], target_genes_list):
            if tarnum > 0:
                other_genes += ge
        other_genes = [gene.name2 for gene in other_genes]  # get gene symbols
        other_genes = set(other_genes)
        motifs_genes[motif] = other_genes

    return bk_genes, motifs_genes


def motif_targenes_de_fisher_test(manorm_pk_fp, motifscan_fp, de_genes_fp, refgenes_fp, is_nearest=False):
    """
    对所有的motifscan中的motif的靶基因做针对给定差异基因的fisher exact检验看是否富集在差异表达基因中
    @param manorm_pk_fp: MAnorm peak文件
    @param motifscan_fp: motifscan tarnum结果文件
    @param refgenes_fp: refseq gene文件
    @param de_genes_fp: 差异基因列表文件或者文件夹
    """
    bkg_gene_set, motifs_genes = get_motifs_target_genes(manorm_pk_fp, motifscan_fp, refgenes_fp, is_nearest)
    pkf_name = os.path.split(manorm_pk_fp)[-1]

    def get_1de_gene_list(one_de_fp):
        de_gene_set = read_gene_list_file(one_de_fp)
        def_name = os.path.split(one_de_fp)[-1]
        fo = open('_'.join([pkf_name[:-4], def_name[:-4], 'fisher_test.xls']), 'w')
        header = \
            '\t'.join(['motif', 'bkg genes', 'motif targenes', 'de genes', 'overlap genes', 'enrich score', 'p-value'])
        fo.write(header + '\n')
        for motif in motifs_genes.keys():
            print '%s\t' % motif,
            out_string = motif + '\t' + make_2geneset_fisher_exact_test(bkg_gene_set, motifs_genes[motif], de_gene_set)
            fo.write(out_string + '\n')
        fo.close()
    if os.path.isfile(de_genes_fp):
        get_1de_gene_list(de_genes_fp)
    elif os.path.isdir(de_genes_fp):
        for fn in os.listdir(de_genes_fp):
            fp = os.sep.join([de_genes_fp, fn])
            try:
                get_1de_gene_list(fp)
            except:
                continue
    print '\n'


def do_peak_target_genes_overlap_fisher_test(manorm_pk_fp, motifscan_fp, refgene_fp, host_motif):
    """
    某个特定的motif的靶基因是否与其他motif的靶基因有富集
    @param manorm_pk_fp: MAnorm peak 文件
    @param motifscan_fp: motifscan tarnum 文件
    @param refgene_fp: refseq基因文件
    @param host_motif: 想要查看的motif
    """
    bkg_gene_set, motifs_genes = get_motifs_target_genes(manorm_pk_fp, motifscan_fp, refgene_fp)

    def do_1host_motif(hm):
        print 'host motif %s targenes fisher test:' % hm
        if hm not in motifs_genes.keys():
            print '@warning: %s not exist!' % hm
            exit(0)

        host_gene_set = motifs_genes[hm]
        for motif in motifs_genes.keys():
            print '%s\t' % motif,
            make_2geneset_fisher_exact_test(bkg_gene_set, motifs_genes[motif], host_gene_set)

    if isinstance(host_motif, list):
        for hm in host_motif:
            do_1host_motif(hm)
    elif isinstance(host_motif, str):
        do_1host_motif(host_motif)


def do_2motif_peaks_overlap_fisher_test(pk_fp, host_motif, motifscan_fp):
    """
    某个特定的motif的peaks是否与其他motif的peaks有富集
    @param pk_fp: 背景peaks
    @param host_motif: 考察的Motif
    @param motifscan_fp: 背景peaks的motifscan结果
    """
    pks, tarnum_dict, motifs = match_manorm_with_motifscan(pk_fp, motifscan_fp)
    bkg_pks_len = len(pks)

    def do_1host_motif(hm):
        print 'motif %s peaks fisher test:' % hm
        host_tarnum = tarnum_dict[hm]
        host_pk_len = __cal_length_over_zero(host_tarnum)
        # other motif
        for moti in motifs:
            other_tarnum = tarnum_dict[moti]
            other_pks_len = __cal_length_over_zero(other_tarnum)
            overlap_tarnum = [h * o for h, o in zip(host_tarnum, other_tarnum)]
            overlap_pks_len = __cal_length_over_zero(overlap_tarnum)
            table = get_fisher_table(bkg_pks_len, host_pk_len, other_pks_len, overlap_pks_len)
            res = fisher_exact(table, alternative='greater')
            enrich_score = cal_enrich_score(bkg_pks_len, host_pk_len, other_pks_len, overlap_pks_len)
            print '%s\t%d\t%d\t%d\t%d\t%.3f\t%s' % \
                  (moti, bkg_pks_len, host_pk_len, other_pks_len, overlap_pks_len, enrich_score, str(res[1]))

    if isinstance(host_motif, str):
        do_1host_motif(host_motif)
    elif isinstance(host_motif, list):
        for hm in host_motif:
            do_1host_motif(hm)


def __cal_length_over_zero(num_list):
    i = 0
    for num in num_list:
        if num > 0:
            i += 1
    return i


def test_motif_targenes_de_fisher_test():
    os.chdir('/mnt/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75')

    pk = \
        '/mnt/MAmotif/2.Processing/4.Histone_Broad_hg19_MAmotif/' \
        'H3K27ac_Rep1_MAmotif/last_version/H3K27ac_H1hesc_VS_K562/' \
        'H1hesc_H3K27ac_Broad_Rep1_distal_peak_MAvalues.xls'
    motifscan = \
        '/mnt/MAmotif/2.Processing/3.Motifscan_result/' \
        'H3K27ac_H1hesc_Rep1_peak_result'
    de_genes = '/mnt/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/DE_genes/' \
               'K562_H1hesc_DEseq_DE_FC_2_p1e-07_genes.txt'
    if 'distal' in pk:
        motif_targenes_de_fisher_test(pk, motifscan, de_genes, hg19ref, True)
    else:
        motif_targenes_de_fisher_test(pk, motifscan, de_genes, hg19ref)


def test_do_2motif_peaks_overlap_fisher_test():
    pk = \
        '/mnt/MAmotif/1.RAWdata/Histone_Broad_hg19/H3K4me3/2.MAnorm/H3K4me3_H1ESC_vs_Helas3_Braod_Rep1_Top12K/' \
        'H1hESC_H3K4me3_Broad_Rep1_Top12K_peaks_MAvalues.xls'
    motifscan = \
        '/mnt/MAmotif/1.RAWdata/Histone_Broad_hg19/H3K4me3/' \
        '3.MotifScan/motifscan_output_H1hESC_H3K4me3_Broad_Rep1_peaks/peak_result.pkl'
    do_2motif_peaks_overlap_fisher_test(pk, key_motif, motifscan)


def test_do_peak_target_genes_overlap_fisher_test():
    os.chdir('F:\\MAmotif.py\\7. primed_naive_hESC')
    pkf = 'hESC_H3K4me3_primed_Top10K_promoter_peak_MAvalues.xls'
    msf = 'hESC_H3K4me3_primed_Top10k.pkl'
    do_peak_target_genes_overlap_fisher_test(pkf, msf, hg19ref, key_motif)


if __name__ == '__main__':
    from constant import hg19_refgenes_file as hg19ref
    from constant import mm9_refgenes_file as mm9ref

    key_motif = 'ZNF263'

    test_motif_targenes_de_fisher_test()
    # test_do_2motif_peaks_overlap_fisher_test()
    pass