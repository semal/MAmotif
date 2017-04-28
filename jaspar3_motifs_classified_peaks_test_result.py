# coding=utf-8
# 将peaks按照是否有某个motif可以分成两类，如果这个motif是细胞特异性的motif，那么它倾向于被对应的关键TF所直接binding，结果是
# 导致这些peak的M值的绝对值更大，当我们拿有此Motif的peaks的M值组和没有此Motif的peaks的M值组做t-test和ranksum-test的时候，
# 观察到的结果就应该是有显著的差异（p值很小）
from optparse import OptionParser
import os
import numpy
from MAmotif_pkg import *
from jaspar3 import jaspar3_name, read_motifscan_result


def get_jaspar3_target_number_list(pk_file, motifscan_result, motifs=jaspar3_name, neg=False):
    """
    read MAnorm peaks and motifscan result, than match the two result.
    """
    # read MAnorm peaks
    if isinstance(pk_file, str):
        print 'read %s' % pk_file
        pk_list = read_MAnorm_peaks(pk_file, neg)
    else:
        pk_list = pk_file

    # read motifscan result
    if isinstance(motifscan_result, list):
        print 'read %s' % motifscan_result
        motifscan_pks = motifscan_result
    else:
        motifscan_pks = read_motifscan_result(motifscan_result, motifs)

    # matching ...
    print 'match MAnorm peaks with motifscan result...'
    tarnum_list = []  # get target number list of matched pks
    match_num = 0
    pk_num = len(pk_list)
    for pk in pk_list:
        match = False
        for mp in motifscan_pks:
            if pk == mp:
                match_num += 1
                os.write(1, '\r%d/%d' % (match_num, pk_num))
                tarnum_list.append(mp.target_number)
                match = True
                break
        if not match:
            tarnum_list.append(None)
            print 'Not Match!'
            pk.prints()

    pk_list = [pk for pk, tarnum in zip(pk_list, tarnum_list) if tarnum is not None]
    tarnum_list = [tarnum for tarnum in tarnum_list if tarnum is not None]
    print '\nend!'
    return pk_list, tarnum_list


def output_jaspar3_motifs_test_result(
        pk_file_path, motifscan_result_path, motifs=jaspar3_name, negative=False, correction_type='benjamin'):
    """
    match peaks with motifscan result once. then output the test result of all jaspar motifs.
    """
    pk_list, tarnum_list = get_jaspar3_target_number_list(pk_file_path, motifscan_result_path, motifs, negative)
    pk_list = [pk for pk, tn in zip(pk_list, tarnum_list) if tn is not None]
    tarnum_list = [tn for tn in tarnum_list if tn is not None]
    pk_array = np.array(pk_list)

    # get target number list with each element list in peaks order,
    # and target number list corresponding with jaspar3 names list
    tarnum = []
    for i, moti in enumerate(motifs):
        tarnum.append(numpy.array([e[i] for e in tarnum_list]))

    # classify pks and do t-test&ranksum-test for each motif
    lines = []
    t_stat, ttest_pvalue, r_stat, rtest_pvalue = [], [], [], []
    for moti, tarnum in zip(motifs, tarnum):
        classifier = MAnormPeaksClassifier(MAnormPeakSet(), Feature())
        yes = pk_array[numpy.where(tarnum > 0)[0]]
        no = pk_array[numpy.where(tarnum == 0)[0]]
        pkset_yes = MAnormPeakSet()
        if yes.size > 0:
            pkset_yes.set_sequences(yes)
        pkset_no = MAnormPeakSet()
        if no.size > 0:
            pkset_no.set_sequences(no)
        classifier.set_feature(pkset_yes, pkset_no)
        line = '%s\t' % moti
        line += '%d\t%f\t%f\t' % (classifier.feature_yes.size, classifier.feature_yes.mean, classifier.feature_yes.std)
        line += '%d\t%f\t%f\t' % (classifier.feature_no.size, classifier.feature_no.mean, classifier.feature_no.std)
        lines.append(line)
        # print line
        ttest = classifier.ttest_feature_classified_peaks()
        t_stat.append(ttest[0])
        ttest_pvalue.append(ttest[1])
        rtest = classifier.ranksum_feature_classified_peaks()
        # rtest = classifier.kstest_feature_classified_peaks()
        r_stat.append(rtest[0])
        rtest_pvalue.append(rtest[1])

    corrected_ttest_pvalue = correct_pvalues(ttest_pvalue, correction_type)
    corrected_rtest_pvalue = correct_pvalues(rtest_pvalue, correction_type)
    max_pvalue = [max(ctp, crp) for ctp, crp in zip(corrected_ttest_pvalue, corrected_rtest_pvalue)]

    # saving test result
    pk_file_name = os.path.split(pk_file_path)[1]
    test_result = open(pk_file_name[:-4].replace('_MAvalues', '') + '_MAmotif_jaspar_output' + '.xls', 'w')
    if correction_type == 'benjamin':
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By Benjamin correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By Benjamin correction\t' \
            'Maximal P-value\n'
    elif correction_type == 'bonferroni':
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By Bonferroni correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By Bonferroni correction\t' \
            'Maximal P-value\n'
    else:
        header = \
            'Motif Name\t' \
            'Target Number\tAverage of Target M-value\tDeviation of Target M-value\t' \
            'Non-target Number\tAverage of Non-target M-value\tDeviation of Non-target M-value\t' \
            'T-test Statistics\tT-test P-value(right-tail)\tT-test P-value By correction\t' \
            'RanSum-test Statistics\tRankSum-test P-value(right-tail)\tRankSum-test P-value By correction\t' \
            'Maximal P-value\n'
    test_result.write(header)

    idx_sorted = np.array(max_pvalue).argsort().tolist()
    for i in idx_sorted:
        line = lines[i]
        print line
        if ttest_pvalue[i] is not None:
            line += str(t_stat[i]) + '\t' + str(ttest_pvalue[i]) + '\t' + str(corrected_ttest_pvalue[i]) + '\t'
            line += str(r_stat[i]) + '\t' + str(rtest_pvalue[i]) + '\t' + str(corrected_rtest_pvalue[i]) + '\t'
            line += str(max_pvalue[i]) + '\n'
        else:
            line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (None, None, None, None, None, None, None)
            # print line
        test_result.write(line)


def correct_pvalues(pvalues, correction_type='benjamin'):
    sorted_indices = np.array(pvalues).argsort().tolist()
    corrected_pvalues = [None] * len(pvalues)
    if correction_type == 'benjamin':
        for rank, i in enumerate(sorted_indices):
            if pvalues[i] is None:
                continue
            corrected_pvalues[i] = len(pvalues) * pvalues[i] / (rank + 1)
            if corrected_pvalues[i] > 1:
                corrected_pvalues[i] = 1
    elif correction_type == 'bonferroni':
        for rank, i in enumerate(sorted_indices):
            if pvalues[i] is None:
                continue
            corrected_pvalues[i] = len(pvalues) * pvalues[i]
            if corrected_pvalues[i] > 1:
                corrected_pvalues[i] = 1
    else:
        corrected_pvalues = pvalues
    return corrected_pvalues


if __name__ == '__main__':
    # -------------------------------command-------------------------------------------------------------------
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
    print negative
    correction = options.correction
    output_jaspar3_motifs_test_result(pk, motif, jaspar3_name, negative=negative, correction_type=correction)

    # ---- test -----
    # pvalues = [0.015, 0.012, 0.011, 0.022]
    # s_pvalues = [0.011, 0.012, 0.015, 0.022]
    # corrected_pvalues = correct_pvalues(pvalues, correction='bonferroni')
    # s_corrected_pvalues = correct_pvalues(s_pvalues, correction='bonferroni')
    # print corrected_pvalues
    # print s_corrected_pvalues

    # ------test for negative test----------------------
    # pk = os.sep.join(['data', 'K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues.xls'])
    # motifscan = os.sep.join(['data', 'K562_H3K27ac_Rep2.pkl'])
    # output_jaspar3_motifs_test_result(pk, motifscan, jaspar3_name, negative=True, correction_type='benjamin')