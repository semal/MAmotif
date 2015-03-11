# coding=utf-8
# 测试用的数据
import os
from MAmotif_pkg import read_refgenes


program_folder = os.path.split(os.path.abspath(__file__))[0]
data_folder = os.sep.join([program_folder, 'data'])

# 各种常用数据格式示例
hg19_refgenes_file = os.sep.join([data_folder, 'modified_hg19_refgene.txt'])
mm9_refgenes_file = os.sep.join([data_folder, 'modified_mm9_refgene.txt'])
manorm_peak_file = os.sep.join([data_folder, 'H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls'])
motifscan_file = os.sep.join([data_folder, 'H3K27ac_H1hesc_Rep2_peak_result'])
gene_fold_change_file = os.sep.join([data_folder, 'H1hesc_K562_genes_fold_change.txt'])
de_gene_list_file = os.sep.join([data_folder, 'K562_H1hesc_DEseq_DE_FC_2_p1e-07_genes.txt'])


def read_motif_file(motif_fp):
    with open(motif_fp) as fi:
        return [line.strip() for line in fi]
jaspar3_file = os.sep.join([data_folder, 'jaspar3_motif_name_list'])
jaspar3_motifs = read_motif_file(jaspar3_file)

# refseq 的 Transcript id 与 gene symbol 的映射字典
hg19_tid_dict = {gene.name: gene.name2 for gene in read_refgenes(hg19_refgenes_file)}
mm9_tid_dict = {gene.name: gene.name2 for gene in read_refgenes(mm9_refgenes_file)}


pass