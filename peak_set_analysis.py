# coding=utf-8
# peaks（蛋白结合区域）之间的分析是很重要的一环
from scipy.stats import fisher_exact
from MAmotif_pkg import *


def get_overlap_num_in_one_chromosome(peaks1_chrm, peaks2_chrm):
    """
    peaks和seqs都在同一条染色体上，计算peaks有多少与seqs重叠
    """
    i = 0
    seqs_start = np.array([seq.start for seq in peaks2_chrm])
    seqs_end = np.array([seq.end for seq in peaks2_chrm])
    overlap_pks = []
    for pk in peaks1_chrm:
        claus = 1.0 * (pk.end - seqs_start) * (seqs_end - pk.start)
        if claus[claus > 0].size > 0:
            i += 1
            overlap_pks.append(pk)
    return i, overlap_pks


def __cal_overlap_num(peaks1, peaks2):
    """
    计算两组peak的overlap
    """
    overlap_num = 0
    overlap_pks = []
    if not isinstance(peaks1, dict):
        peaks1 = peaks1.group_by_chromosomes()
    if not isinstance(peaks2, dict):
        peaks2 = peaks2.group_by_chromosomes()
    for chrm in peaks1.keys():
        if chrm not in peaks2.keys():
            overlap_num += 0
            continue
        overlap_num_chrm, overlap_pks_chrm = \
            get_overlap_num_in_one_chromosome(peaks1[chrm], peaks2[chrm])
        overlap_num += overlap_num_chrm
        overlap_pks += overlap_pks_chrm
    overlap_peak_set = PeakSet()
    overlap_peak_set.set_sequences(overlap_pks)
    overlap_pks = overlap_peak_set.group_by_chromosomes()
    return overlap_num, overlap_pks


def cal_overlap_num(peaks_fp1, peaks_fp2, log2_m_cut=0.):
    """
    输入两组peaks文件，计算它们之间的overlap
    @param peaks_fp1: peak文件1的路径
    @param peaks_fp2: peak文件2的路径
    @return: 重叠的情况（包括peaks）
    """
    if log2_m_cut != 0:
        pks1 = MAnormPeakSet([p for p in read_peaks(peaks_fp1) if p.mvalue >= log2_m_cut])
        pks2 = MAnormPeakSet([p for p in read_peaks(peaks_fp2) if p.mvalue >= log2_m_cut])
    else:
        pks1, pks2 = MAnormPeakSet(read_peaks(peaks_fp1)), MAnormPeakSet(read_peaks(peaks_fp2))
    pks1_overlap_num, pks1_overlap = __cal_overlap_num(pks1, pks2)
    pks2_overlap_num, pks2_overlap = __cal_overlap_num(pks2, pks1)
    return pks1_overlap_num, pks2_overlap_num


def output_group_peaks_overlap_situation(key_pk_fp, inspect_pk_fdp, log2_m_cut=0.):
    """
    一个peak文件和一组peak文件的overlap情况
    @param key_pk_fp: 关键peak文件的路径
    @param inspect_pk_fdp: 被考察的peaks的文件夹
    """
    key_pks = read_peaks(key_pk_fp)
    if log2_m_cut != 0:
        key_pks = [p for p in key_pks if p.mvalue > log2_m_cut]
    key_fn = os.path.basename(key_pk_fp)
    print '%s(%d)' % (key_fn, len(key_pks)),
    for fn in os.listdir(inspect_pk_fdp):
        # print fn
        inspect_pk_fp = os.sep.join([inspect_pk_fdp, fn])
        key_overlap_num, inspect_overlap_num = cal_overlap_num(key_pk_fp, inspect_pk_fp, log2_m_cut)
        print '\t(%d, %d)' % (key_overlap_num, inspect_overlap_num),
    print '\n',


def output_2group_peaks_overlap_situation(key_pk_fdp, inspect_pk_fdp, log2_m_cut=0.):
    """
    看两个文件夹内的peak文件之间的重叠情况
    """
    fn_list = os.listdir(inspect_pk_fdp)
    print '\t%s' % '\t'.join(fn_list)
    for fn in os.listdir(key_pk_fdp):
        key_pk_fp = os.sep.join([key_pk_fdp, fn])
        output_group_peaks_overlap_situation(key_pk_fp, inspect_pk_fdp, log2_m_cut)


def cal_peaks_enrichment(peak_fp, key_peak_fp, simulate_times=1000):
    """
    将一组peaks当做背景参考，称为参考组。然后再给定一组peaks，称为给定组。
    随机模拟一组和给定组一样的peaks，称为模拟组。看给定组和模拟组与参考组的overlap是否有差异
    @param peak_fp:
    @param key_peak_fp:
    @param simulate_times:
    @return:
    """

    print '参考组文件：%s' % peak_fp
    print '给定组文件：%s' % key_peak_fp

    peaks = read_MAnorm_peaks(peak_fp)
    peakset = MAnormPeakSet()
    peakset.set_sequences(peaks)
    group_peaks = peakset.group_by_chromosomes()

    key_peaks = read_seqs(key_peak_fp)
    key_seqs = Sequences()
    key_seqs.set_sequences(key_peaks)
    group_key_seqs = key_seqs.group_by_chromosomes()

    print '参考组peak的数目：%d' % len(peaks)
    print '给定组peak的数目：%d' % len(key_peaks)

    num_overlap = 0
    for chrm in group_key_seqs.keys():
        if chrm not in group_peaks.keys():
            num_overlap += 0
            continue
        num_overlap += get_overlap_num_in_one_chromosome(group_peaks[chrm], group_key_seqs[chrm])[0]
    print '观测到的overlap个数：%d' % num_overlap  # num_overlap是观测到的结果

    num_simulate_overlap = 0
    over_times = 0  # 模拟结果中出现overlap数目大于观测结果的次数
    for i in range(simulate_times):
        temp = 0  # 记录一次模拟中overlap的序列个数
        os.write(1, '\rsimulate_time: %d / %d' % (i + 1, simulate_times))
        for chrm in group_key_seqs.keys():
            if chrm not in hg19_chrlen.keys():
                continue
            if chrm not in group_peaks.keys():
                continue
            simulated_chrm_seqs = generate_a_chromosome_seqs(group_key_seqs[chrm], chrm)  # 对给定组进行随机模拟
            overlap_num_chrm = get_overlap_num_in_one_chromosome(group_peaks[chrm], simulated_chrm_seqs)[0]
            temp += overlap_num_chrm
            num_simulate_overlap += overlap_num_chrm
        sys.stdout.flush()
        if temp >= num_overlap:
            over_times += 1
    num_simulate_overlap /= simulate_times
    print '\n'
    print '随机模拟overlap的个数：%d' % num_simulate_overlap
    print '随机模拟%d次，出现超过观测值%d的数目：%d\n\n' % (simulate_times, num_overlap, over_times)


def get_all_overlap_num(bk_peaks_fp, motif_peaks_fdp, tf_peaks_fdp):
    """
    将所有的motif组的文件放在一个文件夹，将所有tf组的文件放在另一个文件夹。计算所有motif组和tf组两两之间的overlap数目
    """
    bk_peaks = read_MAnorm_peaks(bk_peaks_fp)
    bk_peaks_set = MAnormPeakSet()
    bk_peaks_set.set_sequences(bk_peaks)
    # group_bk_peaks = bk_peaks_set.group_by_chromosomes()
    print '背景组(%s)peaks的数目是：%d' % (bk_peaks_fp, bk_peaks_set.size)

    motifs_peak_num = []
    for motif_file in os.listdir(motif_peaks_fdp):
        motif_file_path = os.sep.join([motif_peaks_fdp, motif_file])
        motif_peaks = read_MAnorm_peaks(motif_file_path)
        motifs_peak_num.append(len(motif_peaks))
    first_line = \
        '\t'.join(['%s(%d)' % (name[:-4], num) for name, num in
                   zip(os.listdir(motif_peaks_fdp), motifs_peak_num)])
    print '\t%s' % first_line

    for tf_file in os.listdir(tf_peaks_fdp):
        tf_file_path = os.sep.join([tf_peaks_fdp, tf_file])
        tf_peaks = read_bed_peaks(tf_file_path)
        tf_peaks_set = PeakSet()
        tf_peaks_set.set_sequences(tf_peaks)
        # group_tf_peaks = tf_peaks_set.group_by_chromosomes()
        tf_overlap_num, tf_overlap_pks = \
            __cal_overlap_num(bk_peaks_set, tf_peaks_set)  # tf peaks与背景组peaks的重叠数目
        # bk_overlap_num = cal_overlap_num(group_tf_peaks, group_bk_peaks)  # tf peaks与背景组peaks的重叠数目
        # print '与背景的overlap数：%d' % bk_overlap_num
        print '%s:%d/%d' % (tf_file.split('.')[0], tf_overlap_num, tf_peaks_set.size),

        for motif_file in os.listdir(motif_peaks_fdp):
            motif_file_path = os.sep.join([motif_peaks_fdp, motif_file])
            motif_peaks = read_MAnorm_peaks(motif_file_path)
            motif_peaks_set = MAnormPeakSet()
            motif_peaks_set.set_sequences(motif_peaks)
            # group_motif_peaks = motif_peaks_set.group_by_chromosomes()
            motif_overlap_num, motif_overlap_pks = \
                __cal_overlap_num(bk_peaks_set, motif_peaks_set)  # motif peaks与背景组peaks的重叠数目
            overlap_num, overlap_pks = \
                __cal_overlap_num(motif_overlap_pks, tf_overlap_pks)  # motif peaks和tf peaks的重叠数目

            # print 'back overlap number: %d' % bk_overlap_num
            # print 'motif overlap number: %d' % motif_overlap_num
            table = [[overlap_num, tf_overlap_num - overlap_num],
                     [motif_overlap_num - overlap_num,
                      bk_peaks_set.size - tf_overlap_num - motif_overlap_num + overlap_num]]
            # print table

            odd_ratio, pvalue = fisher_exact(table, alternative='greater')
            try:
                enrichment_score = \
                    1. * overlap_num / (tf_overlap_num * motif_overlap_num / bk_peaks_set.size)
                print '\t%.7f,%s,%d' % (enrichment_score, str(pvalue), overlap_num),
                # print '\t%.7f' % enrichment_score,
            except ZeroDivisionError:
                enrichment_score = '1'
                print '\t%s,%s,%d' % (enrichment_score, str(pvalue), overlap_num),
                # print '\t%s' % enrichment_score,
        print '\n',


def test_output_two_group_peaks_overlap_situation():
    key_pk_fp = 'F:\\MAmotif\\7. primed_naive_hESC\\7. comparison_with_Broad_peaks\\MAnorm_peaks\\Broad\\' \
                'H1hESC_H3K4me3_Broad_Rep1_peaks_MAvalues.xls'
    inspect_pk_fdp = 'F:\\MAmotif\\7. primed_naive_hESC\\7. comparison_with_Broad_peaks\\MAnorm_peaks\\primed_vs_naive'
    output_group_peaks_overlap_situation(key_pk_fp, inspect_pk_fdp)


def test_output_2group_peaks_overlap_situation():
    key_pk_fdp = 'F:\\MAmotif\\7. primed_naive_hESC\\7.comparison_with_Broad_peaks\\macs_peaks\\Broad_peaks\\cut_peaks'
    inspect_pk_fdp = \
        'F:\\MAmotif\\7. primed_naive_hESC\\7.comparison_with_Broad_peaks\\macs_peaks\\naive_peaks\\cut_peaks'
    output_2group_peaks_overlap_situation(key_pk_fdp, inspect_pk_fdp)

if __name__ == '__main__':
    # test_output_two_group_peaks_overlap_situation()
    # test_output_2group_peaks_overlap_situation()
    os.chdir('F:\\MAmotif\\4. key_peaks_tfs_enrich')
    get_all_overlap_num('H1hesc_H3K27ac_Broad_Rep2_stringent_distal_peak_MAvalues.xls',
                        'H1hESC_H3K27ac_stringent_distal_motifs_peak',
                        '38_tfs_peak')