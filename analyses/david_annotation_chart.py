# coding=utf-8
# 解析david annotation char的输出
import os
from scipy.stats import fisher_exact


class Term(object):
    def __init__(self):
        self.category = ''
        self.term = ''
        self.count = 0
        self.percentage = 0
        self.pvalue = 0
        self.genes = []
        self.list_total = 0
        self.pop_hits = 0
        self.pop_total = 0
        self.fold_enrichment = 0
        self.bonferroni = 0
        self.benjamini = 0
        self.fdr = 0

    def set_term(self, term):
        self.term = term

    def set_genes(self, gene_list):
        self.genes = gene_list

    def set_benjamini(self, value):
        self.benjamini = value


def get_david_annotation_list(david_annotation_chart_file):
    """
    读取david分析后的annotation chart文件，将所有term的信息保存到一个term的list中。
    """
    annotation_list = []
    with open(david_annotation_chart_file) as fi:
        fi.readline()  # skip header
        for line in fi:
            cnt = line.split('\t')
            a_term = Term()
            a_term.set_term(cnt[1])
            gene_list = [e.strip() for e in cnt[5].split(',')]
            a_term.set_genes(gene_list)
            a_term.set_benjamini(float(cnt[11]))
            annotation_list.append(a_term)
    return annotation_list


def read_motif_peaks_target_genes(motif_peak_target_genes_file):
    genes = []
    with open(motif_peak_target_genes_file) as fi:
        for line in fi:
            genes.append(line.strip().upper())
    genes = set(genes)
    print '@info: get %d genes in %s' % (len(genes), os.path.basename(motif_peak_target_genes_file))
    return genes


def calculate_motif_target_genes_enrich_terms(annotation_list, motif_peak_target_genes_file):
    """
    以所有的做此david分析的gene list作为背景，看每个term里面的gene list与给定的基因list的富集情况。
    """
    all_genes = []
    for tm in annotated_list:
        all_genes += tm.genes
    background_genes = set(all_genes)
    print '@info: %d background genes' % len(background_genes)
    # print background_genes
    motif_target_genes = \
        read_motif_peaks_target_genes(motif_peak_target_genes_file).intersection(background_genes)
    print \
        '\t'.join(['term',
                   'overlap',
                   'motif_peak_target_genes',
                   'term_genes',
                   'background',
                   'enrich_score',
                   'pvalue'])
    for tm in annotated_list:
        term_genes = set(tm.genes)
        overlap_genes = term_genes.intersection(motif_target_genes)
        m11, m12 = len(overlap_genes), len(motif_target_genes) - len(overlap_genes)
        m21 = len(term_genes) - len(overlap_genes)
        m22 = len(background_genes) - len(term_genes) - len(motif_target_genes) + len(overlap_genes)
        table = [[m11, m12], [m21, m22]]
        odd_ratio, pvalue = fisher_exact(table, alternative='greater')
        enrich_score = \
            1.0 * len(overlap_genes) / (1.0 * len(motif_target_genes) * len(term_genes) / len(background_genes))
        if pvalue <= 1:
            print '%s\t%d\t%d\t%d\t%d\t%s\t%s' % \
                  (tm.term, len(overlap_genes), len(motif_target_genes), len(term_genes),
                   len(background_genes), str(enrich_score), str(pvalue))


if __name__ == '__main__':
    os.chdir('/mnt/MAmotif/3.Analysis/3.ZNF263TargetGenesAnnotation')
    annotated_list = \
        get_david_annotation_list(
            'Ese14_H3K27ac_LICR_Rep2_promoter_peaks_MAvalues_targenes_GOTERM_MF_result.txt')

    annotated_list = [tm for tm in annotated_list if tm.benjamini < 0.001]
    motif_peak_target_file = \
        'Ese14_H3K27ac_LICR_Rep2_promoter_peaks_MAvalues_exist_MZF1_5-13_targenes.txt'

    calculate_motif_target_genes_enrich_terms(annotated_list, motif_peak_target_file)