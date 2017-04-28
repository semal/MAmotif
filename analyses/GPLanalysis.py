# coding=utf-8
import os
import numpy as np
from numpy.ma import log2


class GPLArray(object):  # 这个类用于处理基因芯片的表达数据
    """
    表达数据的格式如下：
    ID_REF	VALUE	RENORMALIZED_VALUE
    11715100_at	49.952	49.388
    11715101_s_at	133.01	128.18

    不过不同的平台的数据格式是不一样的
    """
    def __init__(self, i, v):
        self.gpl_id = i
        self.value = v

    @staticmethod
    def read_sample_table(sample_file):
        fi = open(sample_file)
        fi.readline()  # skip the header line
        sample_records = \
            [GPLArray(line.split('\t')[0].strip(), float(line.split('\t')[1])) for line in fi]
        return sample_records


class IDMap(object):  # 用于将某基因芯片平台的探针id和refseq transcript id对应起来
    """
    去掉不相关数据后的数据格式如下：
    ID	RefSeq Transcript ID
    11715100_at	NM_003534
    11715101_s_at	NM_003534
    11715102_x_at	NM_003534
    """
    gpl_id = ''  # 对应ID
    ref_id_list = ''  # 对应RefSeq Transcript ID

    def __init__(self, gpl_id, refseq_id):
        self.gpl_id = gpl_id
        self.ref_id_list = refseq_id

    @staticmethod
    def read_map_file(fp):
        res = []
        with open(fp) as fi:
            for line in fi:
                if line.startswith('ID') or line.startswith('#'):
                    print line
                    continue
                cnt = line.strip().split('\t')
                gpl_id = cnt[0].strip()
                ref_id_list = []
                try:
                    cnt[10]
                except IndexError:
                    continue
                if cnt[10] == '---':
                    continue
                for e in cnt[10].strip().split('///'):
                    g = e.strip().split('//')
                    try:
                        g[1]
                    except IndexError:
                        continue
                    if g[1].strip() == 'RefSeq':
                        ref_id_list.append(g[0])
                res.append(IDMap(gpl_id, ref_id_list))
                # print '\t'.join([gpl_id, str(ref_id_list)])
            res_dict = {r.gpl_id: r.ref_id_list for r in res}
            return res_dict


def get_refgenes_expression_value(sample_fp, gpl_map_id_fp):
    new_sample = {}  # key是refseq的id，value是表达值
    sample_records = GPLArray.read_sample_table(sample_fp)
    id_map_dict = IDMap.read_map_file(gpl_map_id_fp)
    for rec in sample_records:
        try:
            for refseq_id in id_map_dict[rec.gpl_id]:
                try:
                    new_sample[refseq_id].append(rec.value)
                except KeyError:
                    new_sample[refseq_id] = []
                    new_sample[refseq_id].append(rec.value)
        except KeyError:
            print '这个探针%s在id映射表里面没有找到！' % rec.gpl_id

    fw = open(sample_fp.replace('.txt', '_refseq.txt'), 'w')
    for key in sorted(new_sample.keys()):
        value = 2 ** ((log2(np.array(new_sample[key]))).mean())
        line = '%s\t%.3f\n' % (key, value)
        fw.write(line)
        # print np.array(new_sample[key]).mean()
        # print new_sample[key]
        # raw_input('暂停一下')


def get_de_genes(gene_de_file, output_fn, log2fc=1.):
    """
    gene expression file由两列组成，第一列是RefSeq_id，第二列是log2FC
    """
    fo = open(output_fn, 'w')
    with open(gene_de_file) as fi:
        fi.readline()  # skip header
        for line in fi:
            cnt = line.strip().split('\t')
            if log2fc > 0:
                if float(cnt[1]) > log2fc:
                    fo.write(cnt[0] + '\n')
            elif log2fc < 0:
                if float(cnt[1]) < log2fc:
                    fo.write(cnt[0] + '\n')
            else:
                print 'The log2(FC) value is wrong!'


def test_get_refgenes_expression_value():
    os.chdir('F:\\MAmotif\\7. primed_naive_hESC\\5.Waizmann_gene_expression')
    map_fp = 'GPL6244-24073.txt'
    for sample in os.listdir('sample'):
        sample_fp = os.sep.join(['sample', sample])
        get_refgenes_expression_value(sample_fp, map_fp)


def test_get_de_genes():
    os.chdir('F:\\MAmotif\\7. primed_naive_hESC\\5.Waizmann_gene_expression')
    get_de_genes('RefSeq_id_log2_FC_primed_vs_naive.txt', 'primed_vs_naive_waizmann_fc-2_genes.txt', -1)


if __name__ == '__main__':
    test_get_de_genes()
    pass