# coding=utf-8
"""Paths to bundled example / reference data."""
import os

from MAmotif_pkg import read_refgenes

program_folder = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(program_folder, 'data')

# Common example data paths
hg19_refgenes_file = os.path.join(data_folder, 'modified_hg19_refgene.txt')
mm9_refgenes_file = os.path.join(data_folder, 'modified_mm9_refgene.txt')
manorm_peak_file = os.path.join(
    data_folder, 'H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls')
motifscan_file = os.path.join(data_folder, 'H3K27ac_H1hesc_Rep2_peak_result')
gene_fold_change_file = os.path.join(data_folder, 'H1hesc_K562_genes_fold_change.txt')
de_gene_list_file = os.path.join(
    data_folder, 'K562_H1hesc_DEseq_DE_FC_2_p1e-07_genes.txt')
jaspar3_file = os.path.join(data_folder, 'jaspar3_motif_name_list')


def read_motif_file(motif_fp):
    with open(motif_fp, encoding='utf-8') as fi:
        return [line.strip() for line in fi]


_jaspar3_motifs = None
_hg19_tid_dict = None
_mm9_tid_dict = None


def get_jaspar3_motifs():
    global _jaspar3_motifs
    if _jaspar3_motifs is None:
        _jaspar3_motifs = read_motif_file(jaspar3_file)
    return _jaspar3_motifs


def get_hg19_tid_dict():
    global _hg19_tid_dict
    if _hg19_tid_dict is None:
        _hg19_tid_dict = {gene.name: gene.name2 for gene in read_refgenes(hg19_refgenes_file)}
    return _hg19_tid_dict


def get_mm9_tid_dict():
    global _mm9_tid_dict
    if _mm9_tid_dict is None:
        _mm9_tid_dict = {gene.name: gene.name2 for gene in read_refgenes(mm9_refgenes_file)}
    return _mm9_tid_dict


# Backwards-compatible lazy attributes via module-level __getattr__ (PEP 562)
def __getattr__(name):
    if name == 'jaspar3_motifs':
        return get_jaspar3_motifs()
    if name == 'hg19_tid_dict':
        return get_hg19_tid_dict()
    if name == 'mm9_tid_dict':
        return get_mm9_tid_dict()
    raise AttributeError('module {!r} has no attribute {!r}'.format(__name__, name))
