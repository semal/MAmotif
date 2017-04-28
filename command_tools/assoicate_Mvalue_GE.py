# coding=utf-8
# 将MAnorm的M值和M值对应的靶基因的基因表达变化进行关联，计算皮尔森相关系数
from scipy.stats import pearsonr
from import_MAmotif import GeneSet, read_refgenes, read_MAnorm_peaks
from io_related import read_fold_change_file


def cal_pearson_correlation(pk_fp, refgene_fp, genes_fold_change_fp):
    """
    将MAnorm peak文件的每一个peak的M值与对应的靶基因的基因变化关联，计算其皮尔森相关系数
    @param refgene_fp: refgene file download from ENCODE
    @param genes_fold_change_fp: a read density result of RNA-seq data
    """
    ref_genes = GeneSet(read_refgenes(refgene_fp))
    pks = read_MAnorm_peaks(pk_fp)
    targenes = [ref_genes.find_target_gene(pk) for pk in pks]
    # 是的peak和靶基因一一对应，将没有靶基因的peak去掉
    pks = [pks[i] for i in range(len(targenes)) if len(targenes[i]) > 0]
    targenes = [genes[0] for genes in targenes if len(genes) > 0]

    gene_log2fcs = read_fold_change_file(genes_fold_change_fp)
    m_values = [pk.mvalue for pk in pks]
    # Gene的name2是gene symbol
    log2fcs = [gene_log2fcs[gene.name2] if gene.name2 in gene_log2fcs.keys() else None for gene in targenes]
    # 将m_values与log2fcs一一映射，将log2fcs中为None映射的m value去掉
    m_values = [m_values[i] for i in range(len(pks)) if log2fcs[i] is not None]
    log2fcs = [value for value in log2fcs if value is not None]

    prs = pearsonr(m_values, log2fcs)
    print prs

    # 输出关联结果
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(m_values, log2fcs)
    plt.xlabel('mvalue')
    plt.ylabel('log2 fold change')
    plt.title('pearsonr: %f' % prs[0])
    plt.savefig(pk_fp[:-4] + '.png')


def command():
    from optparse import OptionParser
    opt_parser = OptionParser()
    opt_parser.add_option('-p', dest='pk', help='cell to cell MAnorm peak file path')
    opt_parser.add_option('-r', dest='refgene', help='refgene file')
    opt_parser.add_option('-f', dest='fold_change', help='cell to cell log2 fold change')
    options, args = opt_parser.parse_args()
    pk = options.pk
    refgene = options.refgene
    fold_change = options.fold_change
    cal_pearson_correlation(pk, refgene, fold_change)
    print 'Done!'


def test_cal_pearson_correlation():
    from constant import manorm_peak_file, hg19_refgenes_file, gene_fold_change_file
    pk_fp = manorm_peak_file
    refgene_fp = hg19_refgenes_file
    fold_change_fp = gene_fold_change_file
    cal_pearson_correlation(pk_fp, refgene_fp, fold_change_fp)


if __name__ == '__main__':
    command()
    # test_cal_pearson_correlation()
    pass
