# coding=utf-8
# 这个脚本通过读取DEseq的输出结果，用户定义p值和fold change的大小来挑选出差异表达基因
import os


# DE = Difference expression 差异表达
def get_DEgenes_from_DEseq_result(deseq_file, limit_pvalue, more_or_less='>', value=2):
    """
    @param more_or_less: 大于或者小于某个指定的fold change值，这个值有value指定
    """
    fi = open(deseq_file)
    fi.readline()  # skip first line 第一行是header

    log2_fcs = {}
    for line in fi:
        cnt = line.strip().split()
        try:
            if limit_pvalue == 1.:
                log2_fcs[cnt[1].replace('"', '')] = float(cnt[-3])
            else:
                if float(cnt[-2]) < limit_pvalue:
                    if more_or_less == '>':
                        if float(cnt[-3]) > value:
                            log2_fcs[cnt[1].replace('"', '')] = float(cnt[-3])
                    elif more_or_less == '<':
                        if float(cnt[-3]) < value:
                            log2_fcs[cnt[1].replace('"', '')] = float(cnt[-3])
        except ValueError, e:
            print e
            continue
    fi.close()
    filename = os.path.split(deseq_file)[1]
    print 'get %d genes in %s!' % (len(log2_fcs.keys()), filename)
    # 输出差异表达基因，每一行是一个gene symbol
    file_name = filename.replace('.txt', 'log2FC_%s_p%s_genes.txt' % (str(value), str(limit_pvalue)))
    file_name = file_name.replace('read_count_', '').replace('output', '')
    handle = open(file_name, 'w')
    for gene_symbol in log2_fcs.keys():
        handle.write('%s\n' % gene_symbol)


def test_get_DEgenes_from_DEseq_result():
    os.chdir('F:\\2.MAmotif\\10. HEK293\\5.Hek293_RNAseq')
    limit_pvalue = 0.01
    up_fc = 1
    down_fc = -1

    get_DEgenes_from_DEseq_result('HEK293_mRNA_UMS_read_count_H1hesc_read_count_DEseq_DE_output.txt',
                                  limit_pvalue, '>', up_fc)
    get_DEgenes_from_DEseq_result('HEK293_mRNA_UMS_read_count_H1hesc_read_count_DEseq_DE_output.txt',
                                  limit_pvalue, '<', down_fc)


if __name__ == '__main__':
    test_get_DEgenes_from_DEseq_result()
