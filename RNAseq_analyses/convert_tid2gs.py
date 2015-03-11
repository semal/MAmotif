# coding=utf-8
import os
from RNA_seq_analysis import RefSeqTools


def __get_tid_gene_symbol_dict(refseq_fp):
    """
    读取refseq基因文件，返回一个transcript id与gene symbol的字典
    """
    rst = RefSeqTools()
    rst.read_refgene_file(refseq_fp)
    return {refgene.name: refgene.name2 for refgene in rst.ref_genes}


def output_gene_symbol(transcript_id_fp, refseq_fp):
    tid2gs_dict = __get_tid_gene_symbol_dict(refseq_fp)
    with open(transcript_id_fp) as fi:
        tid_list = [line.strip() for line in fi]
        gs_list = list(set([tid2gs_dict[tid] for tid in tid_list if tid in tid2gs_dict.keys()]))
        fo = open(transcript_id_fp.replace('.txt', '_GeneSymbol.txt'), 'w')
        [fo.write(gs + '\n') for gs in gs_list]


def test_output_gene_symbol():
    os.chdir('F:\\ChIA-PET\\K562VsMcf7DEgenes')
    for fn in os.listdir('.'):
        if fn.endswith('Tids.txt'):
            output_gene_symbol(fn, hg19_refgenes_file)


if __name__ == '__main__':
    from constant import hg19_refgenes_file

    test_output_gene_symbol()
    # tid_dict = __get_tid_gene_symbol_dict(hg19_refgenes_file)
    # genes = set([tid_dict[t] for t in tid_dict.keys()])
    # print 'gene symbol number: %d' % len(genes)
    pass