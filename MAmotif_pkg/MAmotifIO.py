# coding=utf-8
# 文件读取和输出相关
import os
import pandas as pd
import numpy as np
import sys

from MAmotifPeaks import MAnormPeak, MotifScanPeak, MAnormPeakSet
from MAnormPeaksClassifier import FeaturePvalue, MAnormPeaksClassifier, FeatureMotif, FeaturePromoter, \
    FeatureOtherPeakOverlap, Feature
from sequence import Gene, Peak, Sequence, GeneSet


def coarse(num, keli):
    return int(num / keli) * keli


def read_seqs(file_path):
    """
    seq文件应该最多只有四列，从左到右分别是：chrm, start, end, strand（可以没有）
    @param file_path: seq文件的路径
    @return: 序列列表
    """
    all_seqs = []
    handle = open(file_path)
    for line in handle:
        cnt = line.strip().split()
        chrm = cnt[0].strip()
        start = int(cnt[1].strip())
        end = int(cnt[2].strip())
        try:
            strand = cnt[3].strip()
            seq = Sequence(chrm, start, end, strand)
        except:
            seq = Sequence(chrm, start, end)
        all_seqs.append(seq)
    return all_seqs


def read_peaks(fp):
    try:
        return read_MACS_peaks(fp)
    except:
        try:
            return read_MAnorm_peaks(fp)
        except:
            return None


def read_MAnorm_peaks(file_path, neg=False):
    """
    read MAnorm peaks set into a MAnormPeakSet
    @param file_path:
    @return: list of MAnorm peaks
    """
    # print 'reading MAnorm peaks of %s, start ...' % file_path.split(os.sep)[-1]
    pk_num = 0
    handle = open(file_path)
    pks = []
    for info in handle:
        if not info.startswith('#'):
            cdt = info.strip().split()
            try:
                pk = MAnormPeak(cdt[0], int(cdt[1]), int(cdt[2]), summit=int(cdt[3]))
            except:
                try:
                    pk = MAnormPeak(cdt[0], int(cdt[1]), int(cdt[2]))
                except:
                    continue
            # pk.set_mvalue(float(cdt[4]))
            try:
                if neg:
                    pk.set_mvalue(-round(float(cdt[4]), 1))  # Keep a decimal point
                else:
                    pk.set_mvalue(round(float(cdt[4]), 1))
            except:
                pass
            # pk.set_mvalue(coarse(float(cdt[4]), 0.5))
            try:
                pk.set_avalue(float(cdt[5]))
            except:
                pass
            try:
                pk.set_pvalue(float(cdt[6]))
            except:
                pass
            pks.append(pk)
            pk_num += 1
    handle.close()
    # print 'end! peak number: %d' % pk_num
    return pks


def read_MACS_peaks(file_path):
    handle = open(file_path)
    pks = []
    for info in handle:
        if not info.startswith('#'):
            cdt = info.split()
            if cdt[0].lower() == 'chr':
                continue
            pk = Peak(cdt[0], int(cdt[1]), int(cdt[2]), int(cdt[1]) + int(cdt[4]))
            pks.append(pk)
    handle.close()
    return pks


def read_bed_peaks(fp):
    pks = []
    flag = 0  # 默认是没有chr的
    with open(fp) as fi:
        for line in fi:
            cdt = line.split()
            if flag == 1 or line.startswith('chr'):
                flag = 1
            if flag == 0:
                pk = Peak('chr' + cdt[0], int(cdt[1]), int(cdt[2]))
            else:
                pk = Peak(cdt[0], int(cdt[1]), int(cdt[2]))
            # print pk.prints()
            pks.append(pk)
    return pks


def read_motifscan_peaks(motifscan_file_path, motif_name):
    """
    read MotifScan pandas pickle result
    @param motifscan_file_path:
    @param motif_name:
    @return: FeatureMotif
    """
    motifscan_peaks_list = []

    if motifscan_file_path.endswith('.mat'):  # from shao Matlab MotifScan result
        from scipy import io as sio
        mat = sio.loadmat(motifscan_file_path)
        mat = mat['peak']
        index = None  # index of motif name for finding corresponding target number
        for i, name in enumerate(mat['motif_name'][0][0][0]):
            if name.lower() == motif_name.lower():
                index = i
                break
        for i in xrange(len(mat['chr_id'][0][0])):
            chrm = mat['chr_id'][0][0][i][0][0]
            start = mat['bpstart'][0][0][i][0]
            end = mat['bpend'][0][0][i][0]
            tarnum = mat['motif_tarnum'][0][0][i]
            pk = MotifScanPeak(chrm, start, end)
            pk.set_motif_info(motif_name, tarnum[index])
            motifscan_peaks_list.append(pk)
        return motifscan_peaks_list

    motifscan_result = pd.read_pickle(motifscan_file_path)  # from WangJiawei MotifScan result
    for k, v in motifscan_result.iterrows():
        motifscan_peak = MotifScanPeak(v['chr'], int(v['start']), int(v['end']), int(v['summit']))
        motifscan_peak.set_motif_info(motif_name, int(v[motif_name + '.tarnum']))
        motifscan_peaks_list.append(motifscan_peak)
    return motifscan_peaks_list


def read_refgenes(gene_file):
    """
    read ref genes
    # -------------------------------------------------------------------------------------------
    # field	example	SQL type	info	description

    # bin	637	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
    # name	NM_021010	varchar(255)	values	Name of gene (usually transcript_id from GTF)
    # chrom	chr8	varchar(255)	values	Reference sequence chromosome or scaffold
    # strand	-	char(1)	values	+ or - for strand
    # txStart	6912828	int(10) unsigned	range	Transcription start position
    # txEnd	6914259	int(10) unsigned	range	Transcription end position
    # cdsStart	6912952	int(10) unsigned	range	Coding region start
    # cdsEnd	6914219	int(10) unsigned	range	Coding region end
    # exonCount	2	int(10) unsigned	range	Number of exons
    # exonStarts	6912828,6914047,	longblob	 	Exon start positions
    # exonEnds	6913065,6914259,	longblob	 	Exon end positions
    # score	0	int(11)	range	score
    # name2	DEFA5	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
    # cdsStartStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
    # cdsEndStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
    # exonFrames	1,0,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon
    # --------------------------------------------------------------------------------------------------------
    @param gene_file: gene info file, refgene download from ENCODE
    """
    gene_list = []
    handle = open(gene_file)
    for gene_info in handle:
        ctt = gene_info.split()
        gene = Gene(ctt[1], ctt[12], ctt[2], int(ctt[4]), int(ctt[5]), ctt[3])
        gene_list.append(gene)
    handle.close()
    return gene_list


def read_motifscan_result(motifscan_result):
    """
    read motifscan result with all jaspar3 motifs (not specific motif like MAmotif_pkg.MAmotifIO.read_motifscan_result)
    这个方法是一次性读取所有jaspar3的motifs的target number结果，而不是在MAmotif_pkg中只读取指定的一个motif的target_number结
    果，这样可以只做一次匹配（匹配MAnorm peaks和Motifscan的结果）而将所有peaks按照有没有某个motif进行分类。
    @param motifscan_result: motifscan result from Wangjiawei's peak result file
    @return: motifscan_peak_list, motifs
    """
    print 'reading motifscan result: %s, start!' % motifscan_result.split(os.sep)[-1]
    motifscan = pd.read_pickle(motifscan_result)  # from Wang Jiawei MotifScan result
    motifs = get_motif_names(motifscan)
    motifscan_peak_list = []
    i = 0
    row_num = motifscan.shape[0]
    for k, v in motifscan.iterrows():
        i += 1
        os.write(1, '\r%d/%d' % (i, row_num))
        motifscan_peak = MotifScanPeak(v['chr'], int(v['start']), int(v['end']), int(v['summit']))
        target_number_list = [int(v[name + '.tarnum']) for name in motifs]
        motifscan_peak.set_motif_info(motifs, target_number_list)
        # print motifscan_peak.target_number
        motifscan_peak_list.append(motifscan_peak)
        sys.stdout.flush()
    print '\nend!'
    return motifscan_peak_list, motifs


def get_motif_names(motifscan):
    motifs = []
    for c in motifscan.columns:
        if c.endswith('tarnum'):
            motifs.append(c[:-7])
    return motifs


def match_manorm_with_motifscan(pk_file, motifscan_result, neg_mvalue=False):
    """
    read MAnorm peaks and motifscan result, than match the two result.
    """
    # read MAnorm peaks
    manorm_pks = read_MAnorm_peaks(pk_file, neg_mvalue)
    print 'MAnorm peaks num: %d' % len(manorm_pks)

    # read motifscan result
    motifscan_pks, motifs = read_motifscan_result(motifscan_result)

    # matching ...
    print 'match MAnorm peaks with motifscan result...'
    matched_manorm_pks = []
    matched_tarnum_list = []
    match_num = 0
    for pk in manorm_pks:
        match = False
        for mp in motifscan_pks:
            if pk == mp:
                match_num += 1
                if match_num % 10000 == 0:
                    print '@info:%d matched!' % match_num
                matched_manorm_pks.append(pk)
                matched_tarnum_list.append(mp.target_number)
                # print mp.target_number
                match = True
                break
        if not match:
            print '@warning:Not Match!\t',
            pk.prints()
    print '@info:%d matched!' % match_num

    tarnum_dict = {}
    for i, motif in enumerate(motifs):
        tarnum_dict[motif] = np.array([e[i] for e in matched_tarnum_list])
        # print tarnum_dict[motif]
    print '\nend!'

    return matched_manorm_pks, tarnum_dict, motifs


def output_peaks_target_genes(pk_fp, refgene_fp):
    """
    获取输入的MAnorm peak文件中每一个peak的靶基因，一个peak可以有多个靶基因，多个peak也可能有共同的靶基因
    @param pk_fp: MAnorm peak文件路径
    @param refgene_fp: refseq基因文件路径
    @return: 所有peak的靶基因的总和
    """
    # get genes
    ref_genes = GeneSet(read_refgenes(refgene_fp))
    # get pks
    pks = read_MAnorm_peaks(pk_fp)
    # get pks target genes
    tar_genes = []
    i = 0
    for pk in pks:
        i += 1
        tar_genes += ref_genes.find_target_gene(pk)
        # print '\r%d / %d' % (i, len(pks))
        os.write(1, '\r%d / %d' % (i, len(pks)))
        sys.stdout.flush()
    print '\n'
    # save gene list
    fo = open(pk_fp.split(os.sep)[-1].split('.')[0] + '_targenes.txt', 'w')
    gene_names = set(gene.name2 for gene in tar_genes)
    [fo.write(name + '\n') for name in gene_names]
    fo.close()
    return tar_genes


def classify_MAnorm_pks_by_promoter(pk_fp, refgene_fp):
    """
    将一组peaks按照peak是否与基因的promoter区域有overlap(重叠)将其分成两类
    """
    # get peaks
    pk_file = os.path.split(pk_fp)[1]
    pklist = read_MAnorm_peaks(pk_fp)
    pkset = MAnormPeakSet(pklist)
    # get feature
    feature = FeaturePromoter(read_refgenes(refgene_fp))
    # make a classifier
    classifier = MAnormPeaksClassifier(pkset, feature)
    print 'start to classify peaks by promoter zone ...'
    classifier.classify_by_feature()
    promoter_pkset = classifier.feature_yes
    print 'peak number in promoter: %d' % promoter_pkset.size
    distal_pkset = classifier.feature_no
    print 'peak number in distal: %d' % distal_pkset.size

    # saving classified result
    header = '\t'.join(['#chr', 'start', 'end', 'summit', 'MAnorm_Mvalue', 'MAnorm_Avalue', 'pvalue\n'])
    promoter_file_name = pk_file.replace('peak', 'promoter_peak')
    file_promoter = open(promoter_file_name, 'w')
    file_promoter.write(header)
    [file_promoter.write(pk.tostring()) for pk in promoter_pkset]
    file_promoter.close()
    distal_file_name = pk_file.replace('peak', 'distal_peak')
    file_distal = open(distal_file_name, 'w')
    file_distal.write(header)
    [file_distal.write(pk.tostring()) for pk in distal_pkset]
    file_distal.close()

    print 'Done!'
    return promoter_file_name, distal_file_name


def classify_MAnorm_peaks_by_motif(pks_fp, motifscan_fp, motif):
    """
    using motifscan result of one specific motif to classify peaks.
    """
    # read MAnorm peak file
    pks_file_name = os.path.split(pks_fp)[1]
    pklist = read_MAnorm_peaks(pks_fp)
    pkset = MAnormPeakSet(pklist)

    # get feature
    feature = FeatureMotif(read_motifscan_peaks(motifscan_fp, motif))

    # make a classifier
    classifier = MAnormPeaksClassifier(pkset, feature)
    classifier.classify_by_feature()
    print 'there are %d mismatch!' % feature.mismatch
    motif_yes = classifier.feature_yes
    motif_no = classifier.feature_no

    # saving classified result
    motif = motif.replace('::', '__')
    file_yes = open(pks_file_name[:-4] + '_exist_%s' % motif + pks_file_name[-4:], 'w')
    file_yes.write('\t'.join(['#chr', 'start', 'end', 'summit', 'MAnorm-Mvalue', 'MAnorm-Avalue', 'p-value\n']))
    [file_yes.write(pk.tostring()) for pk in motif_yes]
    file_yes.close()

    file_no = open(pks_file_name[:-4] + '_absent_%s' % motif + pks_file_name[-4:], 'w')
    file_no.write('\t'.join(['#chr', 'start', 'end', 'summit', 'MAnorm-Mvalue', 'MAnorm-Avalue', 'p-value\n']))
    [file_no.write(pk.tostring()) for pk in motif_no]
    file_no.close()


def classify_MAnorm_pks_by_pvalue(pk_fp, pvalue=0.001):
    """
    MAnorm的结果给每个M值都对应了一个pvalue用来衡量变化的显著性，p值越小说明变化越显著，利用给定的p值将一组peaks分成两类
    """
    # get peaks
    pk_file_name = os.path.split(pk_fp)[1]
    pklist = read_MAnorm_peaks(pk_fp)
    pkset = MAnormPeakSet(pklist)
    # get feature
    feature = FeaturePvalue(pvalue)
    # make a classifier
    classifier = MAnormPeaksClassifier(pkset, feature)
    print 'start to classify peaks ...'
    classifier.classify_by_feature()
    promoter_pkset = classifier.feature_yes
    distal_pkset = classifier.feature_no
    # saving classified result
    header = '\t'.join(['#chr', 'start', 'end', 'summit', 'MAnorm_Mvalue', 'MAnorm_Avalue', 'pvalue\n'])
    file_promoter = open(pk_file_name.replace('peak', 'p_less0.001_peak'), 'w')
    file_promoter.write(header)
    [file_promoter.write(pk.tostring()) for pk in promoter_pkset]
    file_promoter.close()
    file_distal = open(pk_file_name.replace('peak', 'p_over0.001_peak'), 'w')
    file_distal.write(header)
    [file_distal.write(pk.tostring()) for pk in distal_pkset]

    file_distal.close()
    print 'Done!'


def other_pks_classified_pks_ttest(pk_fp, peaks_folder, correction_type='benjamin'):
    """
    using DNA binding peaks to classify peaks into 2 groups.
    """
    # get peaks
    pkset = MAnormPeakSet(read_MAnorm_peaks(pk_fp))

    lines = []
    t_stat, ttest_pvalue, r_stat, rtest_pvalue = [], [], [], []
    for fi in os.listdir(peaks_folder):
        # get feature
        fi_path = os.sep.join([peaks_folder, fi])
        dna_binding_pks = read_seqs(fi_path)
        feature = FeatureOtherPeakOverlap(dna_binding_pks)
        # make a classifier
        classifier = MAnormPeaksClassifier(pkset, feature)
        classifier.classify_by_feature()
        line = '%s\t' % fi
        line += '%d\t%f\t%f\t' % (classifier.feature_yes.size, classifier.feature_yes.mean, classifier.feature_yes.std)
        line += '%d\t%f\t%f\t' % (classifier.feature_no.size, classifier.feature_no.mean, classifier.feature_no.std)
        lines.append(line)
        print line
        ttest = classifier.ttest_feature_classified_peaks()
        t_stat.append(ttest[0])
        ttest_pvalue.append(ttest[1])
        rtest = classifier.ranksum_feature_classified_peaks()
        r_stat.append(rtest[0])
        rtest_pvalue.append(rtest[1])

    corrected_ttest_pvalue = correct_pvalues(ttest_pvalue, correction_type)
    corrected_rtest_pvalue = correct_pvalues(rtest_pvalue, correction_type)
    max_pvalue = [max(ctp, crp) for ctp, crp in zip(corrected_ttest_pvalue, corrected_rtest_pvalue)]

    # saving test result
    pk_file_name = os.path.split(pk_fp)[1]
    test_result = open(pk_file_name[:-4].replace('_MAvalues', '') + '_MAmotif_output' + '.xls', 'w')
    header = _get_header(correction_type)
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

    print 'Done!'


def _get_header(correction_type):
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
    return header


def motif_classified_pks_ttest(pk_file_path, motifscan_result_path, negative=False, correction_type='benjamin'):
    """
    match peaks with motifscan result once. then output the test result of all jaspar motifs.
    """
    pk_list, tarnum_dict, motifs = match_manorm_with_motifscan(pk_file_path, motifscan_result_path, negative)
    pk_array = np.array(pk_list)

    # classify pks and do t-test&ranksum-test for each motif
    lines = []
    t_stat, ttest_pvalue, r_stat, rtest_pvalue = [], [], [], []
    for moti in motifs:
        classifier = MAnormPeaksClassifier(MAnormPeakSet(), Feature())
        yes = pk_array[np.where(tarnum_dict[moti] > 0)[0]]
        no = pk_array[np.where(tarnum_dict[moti] == 0)[0]]
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
    header = _get_header(correction_type)
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


def test_classify_MAnorm_pks_by_pvalue():
    classify_MAnorm_pks_by_pvalue(pk)


def test_classify_MAnorm_pks_by_promoter():
    classify_MAnorm_pks_by_promoter(pk, ref)


def test_match_manorm_with_motifscan():
    match_manorm_with_motifscan(pk, ms)


def test_motif_classified_pks_ttest():
    motif_classified_pks_ttest(pk, ms)


def test_read_motifscan_result():
    read_motifscan_result(ms)


if __name__ == '__main__':

    from constant import manorm_peak_file as pk
    from constant import hg19_refgenes_file as ref
    from constant import motifscan_file as ms

    import time
    start = time.clock()
    # test_classify_MAnorm_pks_by_promoter()
    # test_match_manorm_with_motifscan()
    test_motif_classified_pks_ttest()
    # test_read_motifscan_result()

    print time.clock() - start
    pass