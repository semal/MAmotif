# coding=utf-8
# 通过同源基因数据库homology将human和mouse的基因进行映射，同一将human的基因id转成mouse的基因id，这样就能
# 讨论同源基因是否在human和mouse中都存在同一组蛋白修饰。
import os

from scipy.stats import fisher_exact

from MAmotif_fisher_test import get_motifs_target_genes


def do_homology_mapping(homology_file):
    map_res = {}  # key: human gene symbol  -> value: mouse gene symbol
    handle = open(homology_file)
    i = 0
    for line in handle:
        cnt = line.strip().split('\t')
        human_transcript_id = cnt[5]
        mouse_transcript_id = cnt[1]
        # print '%s: %s' % (human_transcript_id, mouse_transcript_id)
        i += 1
        map_res[human_transcript_id] = mouse_transcript_id
        # map_res[mouse_transcript_id] = human_transcript_id
    print 'total homology gene pairs: %d' % i
    handle.close()
    # import cPickle as pkl
    # pkl.dump(map_res, open('human_mouse_homology_mapping_result.pkl', 'wb'))
    return map_res


def fisher_exact_test(human_gene_list, mouse_gene_list, motif_human_genes, motif_mouse_genes, homology_file):
    map_res = do_homology_mapping(homology_file)
    # all using mouse transcript id
    human = read_gene_set(human_gene_list)
    human = set([map_res[h] for h in human if h in map_res.keys()])
    mouse = read_gene_set(mouse_gene_list)
    overlap = human.intersection(mouse)
    motif_human = read_gene_set(motif_human_genes)
    motif_human = set([map_res[h] for h in motif_human if h in map_res.keys()])
    print motif_human
    motif_mouse = read_gene_set(motif_mouse_genes)
    print motif_mouse
    motif_overlap = motif_human.intersection(motif_mouse)
    # table
    m11 = len(list(motif_overlap))
    m12 = len(list(overlap)) - len(list(motif_overlap))
    m21 = len(list(motif_human)) - len(list(motif_overlap))
    m22 = len(human) - len(motif_overlap)
    table = [[m11, m12], [m21, m22]]
    print table
    enrichment_score = 1.0 * len(motif_overlap) / (len(overlap) * len(overlap) / len(human))
    odd_ratio, pvalue = fisher_exact(table)
    print '%f\t%f\t%s\n' % (enrichment_score, odd_ratio, str(pvalue))
    return enrichment_score, pvalue


def fisher_exact_test2(bk_human_pk, motifscan_human, refgene_human,
                       bk_mouse_pk, motifscan_mouse, refgene_mouse,
                       motifs, homology_file):
    """
    检验某motif的peaks在人和老鼠的组蛋白修饰区域同时存在的富集情况
    """
    map_res = do_homology_mapping(homology_file)
    human, motifs_human = get_motifs_target_genes(bk_human_pk, motifscan_human, refgene_human)
    print 'mouse gene number: %d' % len(human)
    human = set([map_res[h] for h in human if h in map_res.keys()])
    mouse, motifs_mouse = get_motifs_target_genes(bk_mouse_pk, motifscan_mouse, refgene_mouse)
    overlap = human.intersection(mouse)
    print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % \
          ('motif name', 'human', 'human_mouse_overlap', 'motif_human', 'motif_human_mouse_overlap',
           'enrichment score', 'odd ratio', 'p-value')
    for motif in motifs:
        motif_human = set([map_res[h] for h in motifs_human[motif] if h in map_res.keys()])
        # motif_mouse = motifs_mouse[motif]
        motif_overlap = motif_human.intersection(overlap)
        # table
        m11 = len(motif_overlap)
        m12 = len(overlap) - len(motif_overlap)
        m21 = len(motif_human) - len(motif_overlap)
        m22 = len(human) - (len(overlap) + len(motif_human) - len(motif_overlap))
        table = [[m11, m12], [m21, m22]]
        # print table
        if len(motif_human) == 0:
            enrichment_score = -1  # 因为motif的peak没有靶基因存在无法计算enrichment score，用-1表示
        else:
            enrichment_score = 1.0 * len(motif_overlap) / (1.0 * len(motif_human) * len(overlap) / len(human))
        odd_ratio, pvalue = fisher_exact(table, alternative='greater')
        print '%s\t%d\t%d\t%d\t%d\t%s\t%f\t%s' % \
              (motif, len(human), len(overlap), len(motif_human), len(motif_overlap),
               str(enrichment_score), odd_ratio, str(pvalue))
        pass


def read_gene_set(gene_list_file):
    genes = []
    handle = open(gene_list_file)
    for g in handle:
        genes.append(g.strip())
    handle.close()
    return set(genes)


def get_read_count_info(read_count_file):
    """
    将read count文件中的信息读入到字典中，字典的keys是基因名
    """
    fi = open(read_count_file)
    header = fi.readline()
    info = {}
    for line in fi:
        gene_name = line.strip().split()[0]
        info[gene_name] = line.replace(gene_name, '').strip()
    fi.close()

    return info, header


def trim_read_count_files(human_read_count_file, mouse_read_count_file, homology_file):
    """
    将人和小鼠的read count文件变得可比，然后作为DEseq的输入
    """
    human, human_header = get_read_count_info(human_read_count_file)
    mouse, mouse_header = get_read_count_info(mouse_read_count_file)
    map_res = do_homology_mapping(homology_file)
    homology_human_genes = map_res.keys()
    homology_mouse_genes = map_res.values()
    # 去除那些不存于同源基因数据库的基因
    [human.pop(key) for key in human.keys() if key not in homology_human_genes]
    [mouse.pop(key) for key in mouse.keys() if key not in homology_mouse_genes]
    # 将人的基因名替换成老鼠的基因名
    for key in human.keys():
        human[map_res[key]] = human.pop(key)
    # 找人和老鼠之间都出现的那部分基因
    human_genes = set(human.keys())
    mouse_genes = set(mouse.keys())
    overlap_genes = human_genes.intersection(mouse_genes)
    # 将重复出现的基因输出
    new_human_read_count_file = open(human_read_count_file.replace('.', '_new.'), 'w')
    new_human_read_count_file.write(human_header)
    new_mouse_read_count_file = open(mouse_read_count_file.replace('.', '_new.'), 'w')
    new_mouse_read_count_file.write(mouse_header)
    for key in overlap_genes:
        human_line = '\t'.join([key, human[key]]) + '\n'
        new_human_read_count_file.write(human_line)
        mouse_line = '\t'.join([key, mouse[key]]) + '\n'
        new_mouse_read_count_file.write(mouse_line)
    new_human_read_count_file.close()
    new_mouse_read_count_file.close()


def read_mouse_gene_symbol(homology_file):
    mouse_gene_symbols = {}
    with open(homology_file) as fi:
        for line in fi:
            cnt = line.strip().split('\t')
            mouse_gene_symbols[cnt[1]] = cnt[3]
    return mouse_gene_symbols


def read_human_gene_symbol(homology_file):
    human_gene_symbols = {}
    with open(homology_file) as fi:
        for line in fi:
            cnt = line.strip().split('\t')
            human_gene_symbols[cnt[5]] = cnt[7]
    return human_gene_symbols


def replace_transcript_id_with_gene_symbol(tid_file, homology_file):
    """
    在DEseq的结果中，用gene symbol替换transcript id
    """
    mouse = read_mouse_gene_symbol(homology_file)
    # human = read_human_gene_symbol(homology_file)
    tfi = open(tid_file)
    header = tfi.readline()
    fi = open(tid_file.replace('.', '_gene_symbol.'), 'w')
    fi.write(header)
    for line in tfi:
        tid = line.split()[1]
        new_line = line.replace(tid, mouse[tid])
        fi.write(new_line)
    tfi.close()
    fi.close()


def test():
    from constant import jaspar3_motifs
    motifs = jaspar3_motifs
    # ------------test1---------------------------------------------------
    homology_file = os.sep.join(['data', 'human_mouse_homolog_gene.txt'])
    # human_gene_list = os.sep.join(['data', 'H1hesc_H3K27ac_Broad_Rep1_promoter_peak_target_genes.txt'])
    # mouse_gene_list = os.sep.join(['data', 'Esb4_H3K27ac_LICR_Rep1_promoter_peak_target_genes.txt'])
    # motif_human_genes = \
    #     os.sep.join(['data', 'H1hesc_H3K27ac_Broad_Rep1_promoter_peak_MAvalues_exist_ZNF263_target_genes.txt'])
    # motif_mouse_genes = \
    #     os.sep.join(['data', 'Esb4_H3K27ac_LICR_Rep1_promoter_peak_MAvalues_exist_ZNF263_target_genes.txt'])
    # fisher_exact_test(human_gene_list, mouse_gene_list, motif_human_genes, motif_mouse_genes, homology_file)

    # -------------test2------------------------------------------------
    bk_human_pk = os.sep.join(['data', 'H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls'])
    motifscan_human = os.sep.join(['data', 'H3K27ac_H1hesc_Rep2_peak_result'])
    refgene_human = os.sep.join(['data', 'modified_hg19_refgene.txt'])
    bk_mouse_pk = os.sep.join(['data', 'Esb4_H3K27ac_LICR_Rep2_promoter_peak_MAvalues.xls'])
    motifscan_mouse = os.sep.join(['data', 'H3K27ac_Esb4_Rep2_peak_result'])
    refgene_mouse = os.sep.join(['data', 'modified_mm9_refgene.txt'])
    fisher_exact_test2(bk_mouse_pk, motifscan_mouse, refgene_mouse,
                       bk_human_pk, motifscan_human, refgene_human,
                       motifs, homology_file)

    # ----------------test3-------------------------------
    # human_read_count_fi = os.sep.join(['data', 'H1hesc_read_count.txt'])
    # mouse_read_count_fi = os.sep.join(['data', 'Esb4_read_count.txt'])
    # trim_read_count_files(human_read_count_fi, mouse_read_count_fi, homology_file)

    # ------------test4-----------------------------------
    # tid_file = os.sep.join(['data', 'Esb4_H1hesc_DEseq_DE_output.txt'])
    # replace_transcript_id_with_gene_symbol(tid_file, homology_file)


if __name__ == '__main__':
    test()
    pass