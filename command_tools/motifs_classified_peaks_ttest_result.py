# coding=utf-8
# 将peaks按照是否有某个motif可以分成两类，如果这个motif是细胞特异性的motif，
# 那么它倾向于被对应的关键TF所直接binding，
# 结果是导致这些peak的M值的绝对值更大，当我们拿有此Motif的peaks的M值组和没有此Motif的peaks的M值组
# 做t-test和rank-sum-test的时候，观察到的结果就应该是有显著的差异（p值很小）
from import_MAmotif import motif_classified_pks_ttest


def command():
    # -------------------------------command-------------------------------------------------------------------
    from optparse import OptionParser

    opt_parser = OptionParser()
    opt_parser.add_option('-p', dest='pk', help='peak file path')
    opt_parser.add_option('-M', dest='motif', help='motifscan result file path')
    opt_parser.add_option('-n', dest='negative', action='store_true', default=False,
                          help='Using negative test of this pk')
    opt_parser.add_option('-c', dest='correction', default='benjamin',
                          help='correction type of pvalues, benjamin or bonferroni, default is no correction')
    options, args = opt_parser.parse_args()
    pk = options.pk
    motif = options.motif
    negative = options.negative
    # print negative
    correction = options.correction
    motif_classified_pks_ttest(pk, motif, negative=negative, correction_type=correction)


def test():
    # ---- test -----
    # pvalues = [0.015, 0.012, 0.011, 0.022]
    # s_pvalues = [0.011, 0.012, 0.015, 0.022]
    # corrected_pvalues = correct_pvalues(pvalues, correction='bonferroni')
    # s_corrected_pvalues = correct_pvalues(s_pvalues, correction='bonferroni')
    # print corrected_pvalues
    # print s_corrected_pvalues

    # ------test for negative test----------------------
    import os
    os.chdir('F:\\MAmotif.py\\7. K562_peaks')
    pk = 'K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues.xls'
    motifscan = 'K562_H3K27ac_Rep2.pkl'
    motif_classified_pks_ttest(pk, motifscan, negative=True, correction_type='benjamin')


if __name__ == '__main__':
    command()