# coding=utf-8
# 当在一个文件夹下有一种细胞对于其他许多细胞的比较结果（MAmotif结果），
# 我们想看某个motif的检验结果中有显著差异的结果数占所有
# 比较数的多少，这个值我们定义为specific_ratio，这个值越大说明这个motif越可能是这个细胞特异的。
import os
import numpy as np


def get_significant_motifs(test_result_file, min_pvalue=1.):
    """
    read test result file
    @param min_pvalue: min pvalue for filtering test results
    @param test_result_file: file path of test result including peaks Mvalue and peaks target genes fold change test
    @return: motif names which are significant
    """
    significant_motifs = []
    with open(test_result_file) as fi:
        for line in fi:
            bonferroni_pvalue = line.strip().split()[13]
            try:
                if float(bonferroni_pvalue) < min_pvalue:
                    significant_motifs.append(line.strip().split()[0])
            except ValueError, e:
                print e
                continue
    return significant_motifs


def stats_test_result(test_result_folder):
    cutted_motifs = []
    for file in os.listdir(test_result_folder):
        fi_path = os.sep.join([test_result_folder, file])
        cutted_motifs += get_significant_motifs(fi_path, min_pvalue=0.001)
    cutted_motifs = np.array(cutted_motifs)
    for motif in list(set(cutted_motifs)):
        print '%s: ' % motif,
        print '%d' % cutted_motifs[cutted_motifs == motif].size


def get_overlap_motifs(mvalue_test_file, fold_change_test_file):
    """
    get intersections of peaks Mvalue significant motifs and peaks target genes fold change significant motifs
    @param mvalue_test_file: file path of peaks mvalue test result
    @param fold_change_test_file: file path of peaks target genes fold change test result
    @return: intersections of two motifs list
    """
    mvalue_significant_motifs = get_significant_motifs(mvalue_test_file, 0.001)
    fold_change_significant_motifs = get_significant_motifs(fold_change_test_file, 0.001)
    return list(set(mvalue_significant_motifs).intersection(set(fold_change_significant_motifs)))


def show_one_motif_test_result(folder_of_comparisons, test_file_name, motif_name, output_file_name):
    """
    @param test_file_name:
    test file name is a file in each comparisons.
    all comparison folder have the same name of this file
    """
    res = open(output_file_name, 'w')
    i = 1
    for comp_f in os.listdir(folder_of_comparisons):
        order = 0  # 此motif的结果在所有motifs的比较结果中的排序
        comp_f_path = os.sep.join([folder_of_comparisons, comp_f])
        if os.path.isdir(comp_f_path):
            comparison_name = comp_f
            try:
                file_path = os.sep.join([folder_of_comparisons, comp_f, test_file_name])
                hd = open(file_path)
            except IOError, e:
                print e
                continue
            for rec in hd:
                if i == 1:
                    header = rec
                    res.write(header.replace('\n', '\torder\n'))
                    i += 1
                if motif_name == rec.split()[0]:
                    res.write(rec.replace(motif_name, comparison_name).replace('\n', '\t%d\n' % order))
                    break
                order += 1
            hd.close()
    res.close()


def show_one_motif_fisher_test_result(folder_of_comparisons, test_file_name,
                                      motif_name, output_file_name):
    """
    @param test_file_name:
    test file name is a file in each comparisons.
    all comparison folder have the same name of this file
    """
    res = open(output_file_name, 'w')
    i = 1
    for comp_f in os.listdir(folder_of_comparisons):
        comp_f_path = os.sep.join([folder_of_comparisons, comp_f])
        if os.path.isdir(comp_f_path):
            comparison_name = comp_f
            try:
                file_path = os.sep.join([folder_of_comparisons, comp_f, test_file_name])
                hd = open(file_path)
            except IOError, e:
                print e
                continue
            for rec in hd:
                if i == 1:
                    header = rec
                    res.write(header)
                    i += 1
                if motif_name == rec.split()[0]:
                    res_rec = '\t'.join([comparison_name, rec])
                    res.write(res_rec)
            hd.close()
    res.close()


def read_motifs_name(motifs_list_file):
    with open(motifs_list_file) as fi:
        return [line.strip() for line in fi]


def extract_one_motif_result_from_MAmotif2(fd, motifs):

    def isexist(motifs, line):
        for motif in motifs:
            if motif in line:
                return True
        return False
    os.chdir(fd)
    for fn in os.listdir('.'):
        with open(fn) as fi:
            for line in fi:
                if isexist(motifs, line):
                    new_line = '%s:' % fn + line
                    print new_line,


# ----------------------------------------- test --------------------------------------------------
def Broad(his, rep, jaspar3_name):
    for motif_name in jaspar3_name:
        # folder_of_comparisons = '/mnt/MAmotif.py/2.Processing/4.Histone_Broad_hg19_MAmotif/%s_%s_MAmotif' % (his, rep)
        folder_of_comparisons = \
            '/mnt/MAmotif.py/2.Processing/4.Histone_Broad_hg19_MAmotif/%s_%s_MAmotif/coarse_graining_0.1' % (his, rep)
        test_file_name = 'H1hesc_%s_Broad_%s_promoter_peak_MAmotif_jaspar_output.xls' % (his, rep)
        output_file_name = 'H1hesc_%s_Broad_%s_%s_promoter_MAmotif_result.xls' % (his, rep, motif_name)
        show_one_motif_test_result(folder_of_comparisons, test_file_name, motif_name, output_file_name)

        # test_file_name = 'H1hesc_%s_Broad_%s_distal_peak_MAmotif_jaspar_output.xls' % (his, rep)
        # output_file_name = 'H1hesc_%s_Broad_%s_%s_distal_MAmotif_result.xls' % (his, rep, motif_name)
        # show_one_motif_test_result(folder_of_comparisons, test_file_name, motif_name, output_file_name)


def LICR(his, rep, jaspar3_name):
    for motif_name in jaspar3_name:
        folder_of_comparisons = \
            '/mnt/MAmotif/2.Processing/4.Histone_wt_mm9_mESC_MAmotif' \
            '/%s_%s_MAmotif' % (his, rep)
        folder_of_comparisons = \
            '/mnt/MAmotif/2.Processing/4.Esb4_and_Ese14_H3K4me3_MAmotif/4.MAmotif_Shao'
        test_file_name = \
            'Ese14_%s_LICR_%s_Top25K_promoter_peaks_MAmotif_jaspar_output.xls' % (his, rep)
        output_file_name = \
            'Ese14_%s_LICR_%s_%s_Top25K_promoter_MAmotif_result.xls' % (his, rep, motif_name)
        test_file_name = \
            'Ese14_H3K4me3_LICR_Rep1_P100_promoter_peaks_MAmotif_jaspar_output.xls'
        output_file_name = \
            'Ese14_H3K4me3_LICR_Rep1_P100_%s_promoter_peaks_MAmotif_result.xls' % motif_name
        show_one_motif_test_result(
            folder_of_comparisons, test_file_name, motif_name, output_file_name)

        test_file_name = \
            'Ese14_%s_LICR_%s_Top25K_distal_peaks_MAmotif_jaspar_output.xls' % (his, rep)
        output_file_name = \
            'Ese14_%s_LICR_%s_%s_Top25K_distal_MAmotif_result.xls' % (his, rep, motif_name)
        test_file_name = \
            'Ese14_H3K4me3_LICR_Rep1_P100_distal_peaks_MAmotif_jaspar_output.xls'
        output_file_name = \
            'Ese14_H3K4me3_LICR_Rep1_P100_%s_distal_peaks_MAmotif_result.xls' % motif_name
        show_one_motif_test_result(
            folder_of_comparisons, test_file_name, motif_name, output_file_name)


def Epi(his, rep, jaspar3_name):
    for motif_name in jaspar3_name:
        folder_of_comparisons = \
            '/mnt/MAmotif.py/2.Processing/4.Histone_Epi_hg19_MAmotif/%s_%s_MAmotif_Epi' % (his, rep)
        if his == 'H3K27ac':
            test_file_name = 'H1hesc_H3K27ac_SAK270_promoter_peak_MAmotif_jaspar_output.xls'
        else:
            test_file_name = 'H1hesc_H3K9ac_SAK68_promoter_peak_MAmotif_jaspar_output.xls'
        output_file_name = \
            'H1hesc_%s_Epi_%s_%s_promoter_MAmotif_result.xls' % (his, rep, motif_name)
        show_one_motif_test_result(
            folder_of_comparisons, test_file_name, motif_name, output_file_name)

        if his == 'H3K27ac':
            test_file_name = 'H1hesc_H3K27ac_SAK270_distal_peak_MAmotif_jaspar_output.xls'
        else:
            test_file_name = 'H1hesc_H3K9ac_SAK68_distal_peak_MAmotif_jaspar_output.xls'

        output_file_name = 'H1hesc_%s_Epi_%s_%s_distal_MAmotif_result.xls' % (his, rep, motif_name)
        show_one_motif_test_result(
            folder_of_comparisons, test_file_name, motif_name, output_file_name)


def test_show_one_motif_test_result():
    motifs = ['ZNF263', 'MZF1_5-13', 'Sox2', 'Pou5f1']
    reps = ['Rep1', 'Rep2']
    for rep in reps:
        for motif in motifs:
            fd = '/mnt/MAmotif/2.Processing/4.Histone_LICR_mm9_Ese14_MAmotif'
            fn = 'Esb4_H3K4me3_LICR_%s_promoter_peaks_MAmotif_jaspar_output.xls' % rep
            fo = 'Esb4_H3K4me3_LICR_%s_promoter_peaks_%s_MAmotif_jaspar_output.xls' % (rep, motif)
            show_one_motif_test_result(fd, fn, motif, fo)


def test_extract_one_motif_result_from_MAmotif2():
    fd = '/mnt/MAmotif/2.Processing/2.Histone_LICR_mm9_Ese14_MAnorm2/H3K9ac/jaspar_output'
    motifs = ['ZNF263', 'Sox2', 'Pou5f1', 'MZF1_5-13']
    extract_one_motif_result_from_MAmotif2(fd, motifs)


if __name__ == '__main__':
    # test_show_one_motif_test_result()
    test_extract_one_motif_result_from_MAmotif2()
    pass