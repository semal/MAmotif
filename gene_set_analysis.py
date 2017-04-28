# coding=utf-8
import os
from io_related import read_gene_list_file
from scipy.stats import fisher_exact


def get_fisher_table(bk_len, a_len, b_len, overlap_len):
    return [[overlap_len, a_len - overlap_len], [b_len - overlap_len, bk_len - a_len - b_len + overlap_len]]


def cal_enrich_score(bk_len, a_len, b_len, overlap_len):
    try:
        return 1. * overlap_len / (a_len * b_len / bk_len)
    except:
        return 0


def get_2gene_list_overlap(genes1_fp, genes2_fp):
    """
    输出两个基因列表重复的部分
    """
    gs1 = read_gene_list_file(genes1_fp)
    gs2 = read_gene_list_file(genes2_fp)
    overlap = gs1.intersection(gs2)
    genes1_fn = os.path.basename(genes1_fp)
    genes2_fn = os.path.basename(genes2_fp)
    # fo = open('_'.join([genes1_fn.split('.')[0], genes2_fn.split('.')[0], 'overlap']) + '.txt', 'w')
    # [fo.write(name + '\n') for name in overlap]
    # fo.close()
    print 'gene set 1(%s): %d' % (genes1_fp, len(gs1))
    print 'gene set 2(%s): %d' % (genes2_fp, len(gs2))
    print 'gene overlap number: %d' % len(overlap)
    # return overlap
    return len(gs1), len(gs2), len(overlap)


def make_2geneset_fisher_exact_test(bkg_gs, gs1, gs2, transcript_id_2symbol=False, refgene_file=''):
    """
    两个基因集合，在给定的背景基因集合下做是否富集的fisher精确检验
    @param bkg_gs: background gene set
    @param gs1: gene1 set
    @param gs2: gene2 set
    @param transcript_id_2symbol: 是否要将transcript id转成gene symbol
    @return: enrich score, pvalue
    """
    if transcript_id_2symbol:
        from MAmotif_pkg import read_refgenes
        refgenes = read_refgenes(refgene_file)
        tid_symbol_dict = {gene.name: gene.name2 for gene in refgenes}

        def convert_transcript_id_2symbol(tid_set):
            return set(tid_symbol_dict[tid] for tid in tid_set)
        bkg_gs, gs1, gs2 = [convert_transcript_id_2symbol(tid_set) for tid_set in [bkg_gs, gs1, gs2]]

    # 构建fisher检验的列联表
    # 去掉非背景中的基因
    gs1, gs2 = gs1.intersection(bkg_gs), gs2.intersection(bkg_gs)
    overlap_gs = gs1.intersection(gs2)
    bk_len = len(bkg_gs)
    gs1_len = len(gs1)
    gs2_len = len(gs2)
    ol_len = len(overlap_gs)
    table = get_fisher_table(len(bkg_gs), len(gs1), len(gs2), len(overlap_gs))
    res = fisher_exact(table, alternative='greater')
    enrich_score = cal_enrich_score(bk_len, gs1_len, gs2_len, ol_len)
    # 打印格式：背景基因数目\t基因列表1数目\t基因列表2数目\t基因列表1和基因列表2重叠的数目\tEnrichment score\tP value\n
    out_string = '%d\t%d\t%d\t%d\t%.3f\t%s' % (bk_len, gs1_len, gs2_len, ol_len, enrich_score, str(res[1]))
    # print out_string
    return enrich_score


def fisher_exact_test_wrapper(background_genes, genes1, genes2, transcript_id_2symbol=False, refgene_file=''):
    """
    三个参数是三个基因列表文件，第一列必须是基因名，每列\t分隔，以#开头的行将不被处理。
    检验在背景基因列表下，基因列表1和基因列表2是否趋向于富集。
    @param background_genes: 作为背景的基因列表文件
    @param genes1: 基因列表文件1
    @param genes2: 基因列表文件2
    """
    bkg_gs, gs1, gs2 = [read_gene_list_file(fp) for fp in [background_genes, genes1, genes2]]
    make_2geneset_fisher_exact_test(bkg_gs, gs1, gs2, transcript_id_2symbol, refgene_file)


def test_get_2gene_list_overlap():
    os.chdir('F:\\3.ChIA-PET\\3.EnhancerStudy\\2.EnhancerDefinedGenes')
    gsf1 = 'F:\\3.ChIA-PET\\3.EnhancerStudy\\2.GenesRelated\\' \
           'K562vsMcf7_DEseq_DE_log2FC_1_p0.01_genes_gene_symbol.txt'
    for fn in os.listdir('.'):
        if 'Symbol' not in fn or 'Over1' not in fn:
            continue
        gsf2 = fn
        l1, l2, overlap = get_2gene_list_overlap(gsf1, gsf2)
        table = get_fisher_table(25334, l1, l2, overlap)
        print fisher_exact(table, alternative='greater')


def temp():
    os.chdir('/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version/TS')
    broad_hESC_highly_expression_genes_fd = '/mnt/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/DE_genes'
    primed_highly_expression_genes = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version/5.DE_genes/' \
                                     'H1hESC_3iL_naive_H1hESC_primed_Ng_DEseq_DE_FC_2_p0.01_genes_gene_symbol.txt'
    primed_gene_set = read_gene_list_file(primed_highly_expression_genes)
    naive_highly_expression_genes = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version/5.DE_genes/' \
                                    'H1hESC_3iL_naive_H1hESC_primed_Ng_DEseq_DE_FC_-2_p0.01_genes_gene_symbol.txt'
    naive_gene_set = read_gene_list_file(naive_highly_expression_genes)

    print '\t'.join(['', 'primed overlap(%d)', 'naive overlap(%d)']) % (len(primed_gene_set), len(naive_gene_set))
    for gene_fn in os.listdir(broad_hESC_highly_expression_genes_fd):
        broad_gene_fp = os.sep.join([broad_hESC_highly_expression_genes_fd, gene_fn])
        broad_gene_set = read_gene_list_file(broad_gene_fp)
        primed_overlap_genes = get_2gene_list_overlap(primed_highly_expression_genes, broad_gene_fp)
        naive_overlap_genes = get_2gene_list_overlap(naive_highly_expression_genes, broad_gene_fp)
        print '\t'.join([gene_fn + '(%d)', '%d', '%d']) % \
              (len(broad_gene_set), len(primed_overlap_genes), len(naive_overlap_genes))


if __name__ == '__main__':
    test_get_2gene_list_overlap()
