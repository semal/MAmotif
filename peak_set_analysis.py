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
    res_string = '%s(%d)' % (key_fn, len(key_pks)),
    res_string = res_string[0]
    for fn in os.listdir(inspect_pk_fdp):
        # print fn
        inspect_pk_fp = os.sep.join([inspect_pk_fdp, fn])
        key_overlap_num, inspect_overlap_num = cal_overlap_num(key_pk_fp, inspect_pk_fp, log2_m_cut)
        mean_overlap = (key_overlap_num + inspect_overlap_num) / 2
        string = '\t%d, %.2f' % \
            (mean_overlap, mean_overlap * 1. / min(len(key_pks), len(read_peaks(inspect_pk_fp)))),
        res_string += string[0]
    string = '\n',
    res_string += string[0]
    return res_string


def output_2group_peaks_overlap_situation(key_pk_fdp, inspect_pk_fdp, log2_m_cut=0.):
    """
    看两个文件夹内的peak文件之间的重叠情况
    """
    fn_list = os.listdir(inspect_pk_fdp)
    res_string = '\t%s' % '\t'.join(fn_list)
    res_string += '\n'
    for fn in os.listdir(key_pk_fdp):
        key_pk_fp = os.sep.join([key_pk_fdp, fn])
        res_string += output_group_peaks_overlap_situation(key_pk_fp, inspect_pk_fdp, log2_m_cut)
    print res_string


def cal_peaks_enrichment(peak_fp, key_peak_fp, simulate_times=1000):
    """
    将一组peaks当做背景参考，称为参考组。然后再给定一组peaks，称为给定组。
    随机模拟一组和给定组一样的peaks，称为模拟组。看给定组和模拟组与参考组的overlap是否有差异
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
            # 对给定组进行随机模拟
            simulated_chrm_seqs = generate_a_chromosome_seqs(group_key_seqs[chrm], chrm)
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


def get_all_overlap_enrich_matirx(bk_peaks_fp, motif_peaks_fdp, tf_peaks_fdp):
    """
    将所有的motif组的文件放在一个文件夹，将所有tf组的文件放在另一个文件夹。
    计算所有motif组和tf组两两之间的overlap数目
    """
    bk_peaks = read_MAnorm_peaks(bk_peaks_fp)
    bk_peaks_set = MAnormPeakSet()
    bk_peaks_set.set_sequences(bk_peaks)
    res_string = ''
    res_string += 'Background(%s) peaks count: %d\n' % (bk_peaks_fp, bk_peaks_set.size)

    motifs_peak_num = []
    motif_peaks_list = []
    for motif_file in os.listdir(motif_peaks_fdp):
        motif_file_path = os.sep.join([motif_peaks_fdp, motif_file])
        motif_peaks = read_MAnorm_peaks(motif_file_path)
        motif_peaks_list.append(motif_peaks)
        motifs_peak_num.append(len(motif_peaks))
    first_line = \
        '\t'.join(['%s(%d)' % (name[:-4], num) for name, num in
                   zip(os.listdir(motif_peaks_fdp), motifs_peak_num)])
    res_string += '\t%s\n' % first_line

    for tf_file in os.listdir(tf_peaks_fdp):
        tf_file_path = os.sep.join([tf_peaks_fdp, tf_file])
        tf_peaks = read_MAnorm_peaks(tf_file_path)
        # try:
        #     tf_peaks = read_bed_peaks(tf_file_path)
        # except:
        #     continue
        tf_peaks_set = PeakSet()
        tf_peaks_set.set_sequences(tf_peaks)
        tf_overlap_num, tf_overlap_pks = \
            __cal_overlap_num(bk_peaks_set, tf_peaks_set)  # tf peaks与背景组peaks的重叠数目
        res_string += '%s:%d/%d' % (tf_file.split('.')[0], tf_overlap_num, tf_peaks_set.size)

        for motif_peaks in motif_peaks_list:
            motif_peaks_set = MAnormPeakSet()
            motif_peaks_set.set_sequences(motif_peaks)
            motif_overlap_num, motif_overlap_pks = \
                __cal_overlap_num(bk_peaks_set, motif_peaks_set)  # motif peaks与背景组peaks的重叠数目
            overlap_num, overlap_pks = \
                __cal_overlap_num(motif_overlap_pks, tf_overlap_pks)  # motif peaks和tf peaks的重叠数目

            table = [[overlap_num, tf_overlap_num - overlap_num],
                     [motif_overlap_num - overlap_num,
                      bk_peaks_set.size - tf_overlap_num - motif_overlap_num + overlap_num]]

            odd_ratio, pvalue = fisher_exact(table, alternative='greater')
            try:
                enrichment_score = \
                    1. * overlap_num / (tf_overlap_num * motif_overlap_num / bk_peaks_set.size)
                res_string += '\t%.7f,%s,%d' % (enrichment_score, str(pvalue), overlap_num)
            except ZeroDivisionError:
                enrichment_score = '1'
                res_string += '\t%s,%s,%d' % (enrichment_score, str(pvalue), overlap_num)
        res_string += '\n'
    print res_string


def get_stringent_distal_peaks(pkfp, H3K4me3_fp):
        """
        用K4me3的peaks去除掉K27ac distal区域与K4me3的有重叠的部分
        """
        # 读取并合并k4me3的两个replicate
        k4me3_rep1 = read_MAnorm_peaks(H3K4me3_fp)
        seqs1 = Sequences()
        seqs1.set_sequences(k4me3_rep1)

        # k4me3_rep2 = read_MACS_peaks(os.sep.join(['data', 'H1hesc_H3k4me3_Rep2_macs13_peaks_nomodel.xls']))
        # seqs2 = Sequences()
        # seqs2.set_sequences(k4me3_rep2)

        seqs = seqs1
        group_seqs = seqs.group_by_chromosomes()

        # 读取k27ac distal的peaks
        k27ac_distal = read_MAnorm_peaks(pkfp)
        peaks = MAnormPeakSet()
        peaks.set_sequences(k27ac_distal)
        group_peaks = peaks.group_by_chromosomes()
        # 去除到与k4me3有overlap的peaks
        left_peaks = []
        for chrm in group_peaks.keys():
            if chrm not in group_seqs.keys():
                continue
            seqs_start = np.array([seq.start for seq in group_seqs[chrm]])
            seqs_end = np.array([seq.end for seq in group_seqs[chrm]])
            for pk in group_peaks[chrm]:
                claus = 1.0 * (pk.end - seqs_start) * (seqs_end - pk.start)
                if claus[claus > 0].size == 0:
                    left_peaks.append(pk)

        # 输出剩下的peaks
        left_peak_set = MAnormPeakSet()
        left_peak_set.set_sequences(left_peaks)
        file2save = pkfp.replace('distal', 'stringent_distal')
        fi = open(file2save, 'w')
        [fi.write(pk.tostring()) for pk in left_peak_set]


def get_overlap_flags(bk_pk_fp, tf_pk_fdp, res2save):
    """
    将H3K27ac promoter的peak与一组其他peaks做overlap，
    记录H3K27ac promoter每一个peak的其他peaks的overlap情况
    """
    bk_peaks = read_MAnorm_peaks(bk_pk_fp)
    bk_peaks_set = MAnormPeakSet()
    bk_peaks_set.set_sequences(bk_peaks)
    group_bk_peaks = bk_peaks_set.group_by_chromosomes()
    print 'background peaks count: %d' % bk_peaks_set.size

    tf_overlap_flags = {}
    for tf_file in os.listdir(tf_pk_fdp):
        print 'calculate %s ...' % tf_file
        tf_overlap_flags[tf_file] = []
        tf_file_path = os.sep.join([tf_pk_fdp, tf_file])
        tf_peaks = read_MAnorm_peaks(tf_file_path)
        tf_peaks_set = MAnormPeakSet()
        tf_peaks_set.set_sequences(tf_peaks)
        group_tf_peaks = tf_peaks_set.group_by_chromosomes()

        for chrm in group_bk_peaks.keys():

            if chrm not in group_tf_peaks.keys():
                for pk in group_bk_peaks[chrm]:
                    tf_overlap_flags[tf_file].append(0)
                continue

            tf_starts = np.array([pk.start for pk in group_tf_peaks[chrm]])
            tf_ends = np.array([pk.end for pk in group_tf_peaks[chrm]])
            for pk in group_bk_peaks[chrm]:
                claus = 1.0 * (pk.end - tf_starts) * (tf_ends - pk.start)
                if claus[claus > 0].size > 0:
                    tf_overlap_flags[tf_file].append(1)
                else:
                    tf_overlap_flags[tf_file].append(0)
    fi = open(res2save, 'w')
    header = '\t'.join(['', '', ''] + tf_overlap_flags.keys()) + '\n'
    fi.write(header)
    for i in range(bk_peaks_set.size):
        line = \
            '%s\t%d\t%d\t' % (bk_peaks_set[i].chrm, bk_peaks_set[i].start, bk_peaks_set[i].end) + \
            '\t'.join([str(tf_overlap_flags[key][i]) for key in tf_overlap_flags.keys()]) + '\n'
        fi.write(line)
    fi.close()


def test_output_group_peaks_overlap_situation():
    key_pk_fp = '/mnt/MAmotif/1.RAWdata/3.ValidatedUsingPeaksAndReads/K27ac_wt_FDR5_peaks.xls'
    inspect_pk_fdp = '/mnt/MAmotif/2.Processing/1.Histone_LICR_mm9_peaks/cut_peaks/H3K27ac/Rep1'
    output_group_peaks_overlap_situation(key_pk_fp, inspect_pk_fdp)


def test_output_2group_peaks_overlap_situation():
    os.chdir('/mnt/MAmotif/3.Analysis/2.mm9KeyPeaksTFsEnrich/Ese14')
    key_pk_fdp = 'MotifExistPeaksInStringentDistal'
    inspect_pk_fdp = 'MotifExistPeaksInStringentDistal'
    output_2group_peaks_overlap_situation(key_pk_fdp, inspect_pk_fdp)


def test_get_all_overlap_enrich_matrix():
    os.chdir('/mnt/MAmotif/3.Analysis/2.mm9KeyPeaksTFsEnrich')
    bk_fp = 'Esb4/Esb4_H3K27ac_LICR_Rep2_stringent_distal_peak_MAvalues.xls'
    mfpd = 'Esb4/MotifExistPeakInStringentDistal'
    tffpd = 'Esb4/MotifExistPeakInStringentDistal/'
    get_all_overlap_enrich_matirx(bk_fp, mfpd, tffpd)


def test_get_stringent_distal_peaks():
    pkfp = \
        '/mnt/MAmotif/3.Analysis/2.mm9KeyPeaksTFsEnrich/Ese14/' \
        'Ese14_H3K27ac_LICR_Rep2_distal_peaks_MAvalues.xls'
    H3K4me3_fp = '/mnt/MAmotif/2.Processing/4.Esb4_and_Ese14_H3K4me3_MAmotif/2.MAnorm/' \
                 'H3K4me3_Ese14_vs_Cbellum_Rep1_P100/Ese14_H3K4me3_LICR_Rep1_P100_peaks_MAvalues.xls'
    get_stringent_distal_peaks(pkfp, H3K4me3_fp)


def test_get_overlap_flags():
    os.chdir('/mnt/MAmotif/3.Analysis/2.mm9KeyPeaksTFsEnrich/Ese14')
    pf = 'Ese14_H3K27ac_LICR_Rep2_stringent_distal_peaks_MAvalues.xls'
    tf = '/mnt/MAmotif/3.Analysis/2.mm9KeyPeaksTFsEnrich/Ese14/MotifExistPeaksInStringentDistal/'
    get_overlap_flags(pf, tf,
                      'Ese14_H3K27ac_stringent_distal_4col_overlap_flags_with_other_peaks.xls')


if __name__ == '__main__':
    # test_output_group_peaks_overlap_situation()
    # test_output_2group_peaks_overlap_situation()
    test_get_all_overlap_enrich_matrix()
    # test_get_stringent_distal_peaks()
    # test_get_overlap_flags()