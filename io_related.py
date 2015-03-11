# coding=utf-8
def read_fold_change_file(genes_fold_change_fp):
    """
    读取基因表达变化文件，文件的格式如下:
        H1hesc_mean	K562_mean	fold_change	log2_fold_change
    DEFA5	0.020756	0.009441	2.19857720034	1.1365701919
    EPB41L5	3.059781	1.238889	2.46977762682	1.30438115066
    MCM5	72.070998	38.606702	1.86680018048	0.90056751189
    第一列是基因名，第二列是分子表达量，第三列是分母表达量，第四列是分子表达量与分母表达量的比值，第五列是log2的fold change
    @param genes_fold_change_fp: 基因表达变化文件的路径
    @return: dict, 返回基因文件第一列基因名和第五列log2 fold change作为键值对的字典
    """
    with open(genes_fold_change_fp) as fi:
        fi.readline()  # 跳过文件头行
        return {line.split('\t')[0].strip(): float(line.strip().split()[4]) for line in fi}


def read_gene_list_file(gene_fp):
    """
    第一列必须是基因名，每列\t分隔，以#开头的行将不被处理。
    """
    with open(gene_fp) as fi:
        return set(line.split('\t')[0].strip() for line in fi if not line.startswith('#'))


def output_list(list, file_name):
    """
    将list作为文本输出
    @param list: array like list or set
    @param file_name: 输出文本的文件名
    """
    fo = open(file_name, 'w')
    [fo.write(str(line) + '\n') for line in list]


def test_read_fold_change_file():
    from constant import gene_fold_change_file
    tt = read_fold_change_file(gene_fold_change_file)
    for key in tt.keys():
        print '%s: %f' % (key, tt[key])


if __name__ == '__main__':
    test_read_fold_change_file()