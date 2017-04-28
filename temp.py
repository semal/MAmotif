# coding=utf-8
from numpy.random.mtrand import normal
from scipy.stats import fisher_exact
from scipy.stats.distributions import norm

from MAmotif_pkg import *
import constant
import io_related


def temp1():
        """
        用K4me3的peaks去除掉K27ac distal区域与K4me3的有重叠的部分
        """
        # 读取并合并k4me3的两个replicate
        k4me3_rep1 = read_MAnorm_peaks(os.sep.join(['data', 'K4me3_wt_pcut10_bound_region.txt']))
        seqs1 = Sequences()
        seqs1.set_sequences(k4me3_rep1)

        # k4me3_rep2 = read_MACS_peaks(os.sep.join(['data', 'H1hesc_H3k4me3_Rep2_macs13_peaks_nomodel.xls']))
        # seqs2 = Sequences()
        # seqs2.set_sequences(k4me3_rep2)

        seqs = seqs1
        group_seqs = seqs.group_by_chromosomes()

        # 读取k27ac distal的peaks
        k27ac_distal = read_MAnorm_peaks(os.sep.join(['data', 'Esb4_H3K27ac_LICR_Rep2_distal_peak_MAvalues.xls']))
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
        fi = open(
            'data/mouse_distal_motif_peaks/data/mouse_stringent_distal_motif_peaks/'
            'Esb4_H3K27ac_LICR_Rep2_stringent_distal_peak_MAvalues.xls', 'w')
        [fi.write(pk.tostring()) for pk in left_peak_set]


def temp2():
    """
    将H3K27ac promoter的peak与一组其他peaks做overlap，记录H3K27ac promoter每一个peak的其他peaks的overlap情况
    """
    bk_peaks_file = os.sep.join(['data', 'Esb4_H3K27ac_LICR_Rep2_stringent_distal_peak_MAvalues.xls'])
    bk_peaks = read_MAnorm_peaks(bk_peaks_file)
    bk_peaks_set = MAnormPeakSet()
    bk_peaks_set.set_sequences(bk_peaks)
    group_bk_peaks = bk_peaks_set.group_by_chromosomes()
    print '背景组peaks的数目是：%d' % bk_peaks_set.size

    tf_peaks_folder = os.sep.join(['key_peaks', 'mouse_TF_peaks'])
    tf_overlap_flags = {}
    for tf_file in os.listdir(tf_peaks_folder):
        print 'calculate %s ...' % tf_file
        tf_overlap_flags[tf_file] = []
        tf_file_path = os.sep.join([tf_peaks_folder, tf_file])
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

    fi = open('Esb4_H3K27ac_stringent_distal_overlap_flags_with_other_peaks.xls', 'w')
    header = '\t'.join(['', '', ''] + tf_overlap_flags.keys()) + '\n'
    fi.write(header)
    for i in range(bk_peaks_set.size):
        line = \
            '%s\t%d\t%d\t' % (bk_peaks_set[i].chrm, bk_peaks_set[i].start, bk_peaks_set[i].end) + \
            '\t'.join([str(tf_overlap_flags[key][i]) for key in tf_overlap_flags.keys()]) + '\n'
        fi.write(line)
    fi.close()


def temp3():
    os.chdir('F:\\MAmotif.py\\4. ZNF263_target_genes')
    from analyses.human_mouse_homology_analyses import do_homology_mapping
    map_res = do_homology_mapping('human_mouse_homolog_gene.txt')
    background_genes = open('H1hesc_H3K27ac_Rep2_promoter_peaks_target_genes.txt')
    background_gene_set = set(map_res[gene.strip()] for gene in background_genes if gene.strip() in map_res.keys())
    znf263_genes = open('H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues_exist_ZNF263_target_genes.txt')
    znf263_gene_set = set(map_res[gene.strip()] for gene in znf263_genes if gene.strip() in map_res.keys())
    deseq = open('H1hesc_VS_Esb4_DEseq.txt')
    deseq.readline()
    de_gene_set = \
        set(gene.split('\t')[0] for gene in deseq
            if float(gene.split('\t')[5]) > 2 and float(gene.split('\t')[6]) < 0.001)
    de_gene_set = set(gene for gene in de_gene_set)
    background_overlap = background_gene_set.intersection(de_gene_set)
    znf263_overlap = znf263_gene_set.intersection(de_gene_set)
    m11 = len(znf263_overlap)
    m12 = len(background_overlap) - len(znf263_overlap)
    m21 = len(znf263_gene_set) - len(znf263_overlap)
    m22 = len(background_gene_set) - len(background_overlap) - m21
    table = [[m11, m12], [m21, m22]]
    print table
    print fisher_exact(table, alternative='greater')


def temp4():
    os.chdir('F:\\MAmotif.py\\7. primed_naive_hESC')

    def read_genes(fp):
        with open(fp) as fi:
            genes = set([line.split()[0].strip() for line in fi])
            return genes

    def read_sam_genes(fp):
        with open(fp) as fi:
            fi.readline()
            genes = set([line.split('\t')[2].strip() for line in fi])
            return genes

    background = read_genes('H3K4me3_background_genes.txt')
    # de_genes = read_sam_genes('up_regulated_genes_fc2.0_fdr1.txt')
    de_genes = read_genes('primed_genes_over_logfc_1.txt')
    znf263_genes = read_genes('primed_znf263_absent_targenes.txt')
    overlap = de_genes.intersection(znf263_genes)
    from scipy.stats import fisher_exact
    bk_len = len(background)
    de_len = len(de_genes)
    znf263_len = len(znf263_genes)
    ol_len = len(overlap)
    table = [[ol_len, de_len - ol_len], [znf263_len - ol_len, bk_len - de_len - znf263_len + ol_len]]
    print fisher_exact(table, alternative='greater')
    pass


def temp5():
    os.chdir('F:\\MAmotif.py\\7. primed_naive_hESC')
    fi = open('SAM_OUTPUT_UP_regulated_tabl_ctrl_Naive_samp_Primed.txt')
    header = fi.readline()
    fi = fi.readlines()
    given_fc = 5  # 从sam结果中挑选最小值为given_fc的条目
    given_fdr_list = [1, 2, 5]

    def tt(fi, given_fc, given_fdr):
        fo = open('up_regulated_genes_fc%.1f_fdr%d.txt' % (given_fc, given_fdr), 'w')
        fo.write(header)
        for line in fi:
            cnt = line.split()
            fc = float(cnt[-2])
            fdr = float(cnt[-1])
            if fc > given_fc and fdr < given_fdr:
                fo.write(line)

    [tt(fi, given_fc, given_fdr) for given_fdr in given_fdr_list]


def temp6():
    def read_list(fp):
        with open(fp) as fi:
            return set([line.strip() for line in fi])
    fd = 'F:\\MAmotif\\7. primed_naive_hESC\\top30_motifs_overlap\\de_genes\\down'
    sets = {}
    for fn in os.listdir(fd):
        fp = os.sep.join([fd, fn])
        sets[fn] = read_list(fp)
    keys = sets.keys()
    print '\t'.join([''] + keys)
    common = sets[keys[0]]
    for ko in keys:
        common = common.intersection(sets[ko])
        print '%s(%d)' % (ko, len(sets[ko])),
        for ki in keys:
            print '\t%d' % len(sets[ki].intersection(sets[ko])),
        print '\n',
    print 'common genes size: %d(%s)' % (len(common), str(common))


def temp7():
    fp = 'F:\\MAmotif\\7. primed_naive_hESC\\5.MIT_gene_expression\\DE_genes\\' \
         'renormalized_MIT_WIBR2_hESC_primed_vs_naive_fc.txt'
    with open(fp) as fi:
        header = fi.readline()
        print header
        for line in fi:
            cnt = line.strip().split('\t')
            if float(cnt[0]) > 50 and float(cnt[1]) < -2:
                print '\t'.join([cnt[3], cnt[2]])


def temp8():
    fc_file = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version/5.DE_genes' \
              '/MIT_paper_offered_gene_fc_naive_vs_primed.txt'
    max_pvalue = 0.01
    min_log2_fc = -1
    de_genes = []
    # 上面的gene fold change共3列，分别是Log2_FC_Naive vs. hESM 	Pvalue_Naive vs. hESM	Gene_sy bol
    with open(fc_file) as fi:
        fi.readline()
        for line in fi:
            cnt = line.strip().split('\t')
            log2_fc = float(cnt[0])  # log2 FC naive vs hESC
            pvalue = float(cnt[1])
            if not (pvalue < max_pvalue and log2_fc < min_log2_fc):
                continue
            if '///' not in cnt[2]:
                gene_symbol = cnt[2].strip()
                de_genes.append(gene_symbol)
            else:
                gene_symbols = cnt[2].strip().split('///')
                gene_symbols = [gene.strip() for gene in gene_symbols]
                de_genes += gene_symbols
    de_genes = set(de_genes)
    for gene in de_genes:
        print gene


def temp9():
    os.chdir('/mnt/MAmotif/2.Processing/7.primed_vs_naive/5.DE_genes/')
    pkfp = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/2.MAnorm/H3K4me3_hESC_WIBR2_primed_vs_naive_MIT/' \
           'hESC_WIBR2_naive_H3K4me3_MIT_peaks_MAvalues.xls'

    pklist = read_MAnorm_peaks(pkfp)
    de_pkset = MAnormPeakSet([pk for pk in pklist if pk.mvalue < -1 and pk.pvalue < 0.01])  # 选取差异结合区域

    targenes = []
    refgenes = GeneSet(read_refgenes(constant.hg19_refgenes_file))
    for pk in de_pkset:
        targenes += refgenes.find_target_gene(pk)

    targenes = set([gene.name2 for gene in targenes])
    io_related.output_list(targenes, 'MIT_naive_H3K4me3_M_-1_p0.01_peaks_targenes.txt')


def temp10():
    os.chdir('/mnt/MAmotif/2.Processing/7.primed_vs_naive/5.DE_genes')
    refgenes = GeneSet(read_refgenes(constant.hg19_refgenes_file))
    m1_pkfp = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/2.MAnorm/H3K27ac_H1hESC_primed_vs_naive_Ng/' \
              'H1hESC_primed_H3K27ac_Ng_peaks_MAvalues.xls'
    m1_pks = MAnormPeakSet(read_MAnorm_peaks(m1_pkfp))
    m2_pkfp = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/2.MAnorm/H3K4me3_H1hESC_primed_vs_naive_Ng/' \
              'H1hESC_primed_H3K4me3_Ng_peaks_MAvalues.xls'
    m2_pks = MAnormPeakSet(read_MAnorm_peaks(m2_pkfp))
    de_genes = io_related.read_gene_list_file('Ng_primed_up_regulated_genes_log2FC_1_p0.01_genes_gene_symbol.txt')
    fo = open('Ng_primed_vs_naive_features_of_genes.txt', 'w')
    i = 0
    for gene in refgenes:
        i += 1
        if i % 1000 == 0:
            print i
        m1 = m1_pks.find_target_pks(gene)
        m2 = m2_pks.find_target_pks(gene)
        if m1 == [] or m2 == []:
            continue
        is_de = 1 if gene.name2 in de_genes else 0
        # print m1
        # print m2
        # print is_de
        line = '\t'.join(['%s', '%s', '%s', '%s', '%d\n']) % \
              (gene.name, gene.name2,
               str(max([m.mvalue for m in m1])) if m1 != [] else 'None',
               str(max([m.mvalue for m in m2])) if m2 != [] else 'None', is_de)
        fo.write(line)
        # print line
    pass


def temp11():
    header = 'name\tmotif\tbkg genes\tmotif targenes\tde genes\toverlap genes\tenrich score\tp-value'
    print header
    motifs = ['MZF1_5-13', 'Pou5f1', 'Sox2']
    os.chdir('F:\MAmotif\8. jaspar_motifs_targenes_DE_fisher_test\Esb4_DE_genes_fisher_test')
    for motif in motifs:
        for fn in os.listdir('Esb4_promoter'):
            with open(os.sep.join(['Esb4_promoter', fn])) as fi:
                for line in fi:
                    if motif in line:
                        print '%s\t%s' % (fn, line.strip())


def download_38_tf_peaks():
    first_fp = 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1505nnn'  # 一级目录
    gsm_prefix = 'GSM1505'
    gsm_suffix = range(615, 819, 1)
    for suffix in gsm_suffix:
        second_fd = ''.join([gsm_prefix, str(suffix)])  # 二级目录名
        download_fp = '/'.join([first_fp, second_fd, 'suppl', '*.peak.txt.gz'])  # 下载路径
        cmd = 'wget -c %s' % download_fp  # 下载命令
        print cmd

if __name__ == '__main__':
    # temp11()
    download_38_tf_peaks()
