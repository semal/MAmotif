# coding=utf-8
# 计算refgene中每个基因的read count。所有样本都放在一个文件夹里面作为输入。
from optparse import OptionParser
from RNA_seq_analysis import *


def output_read_count(cell_folder, refgenes_file):
    """
    using multiple process to get rpkm, process number is 4.
    saving the result in a txt file named by cell_folder name
    """
    if cell_folder.endswith(os.sep):
        cell_folder = cell_folder[:-1]

    from joblib import Parallel, delayed
    fipath_list = [os.sep.join([cell_folder, fi]) for fi in os.listdir(cell_folder) if fi.endswith('.bed')]
    Parallel(n_jobs=4)(delayed(run2)(fi, refgenes_file) for fi in fipath_list)

    # [run(fi, refgene_file) for fi in fipath_list]
    new_fipath_list = [fipath[:-4] + '_read_count.txt' for fipath in fipath_list]
    # print new_fipath_list

    result = open(cell_folder + '_read_count.txt', 'w')
    header = '\t'.join([''] + [fi for fi in os.listdir(cell_folder) if fi.endswith('.bed')]) + '\n'
    result.write(header)
    rpkms = [[] for _ in range(len(fipath_list))]
    refgenes = []
    for i, fipath in enumerate(new_fipath_list):
        fi = open(fipath)
        fi.readline()  # skip the header line
        for line in fi:
            rpkms[i].append(line.split()[1])
            if i == 0:
                refgenes.append(line.split()[0])
        fi.close()
    # print rpkms
    for i, refgene in enumerate(refgenes):
        line = '\t'.join([refgene] + [rpkm[i] for rpkm in rpkms]) + '\n'
        result.write(line)
    result.close()
    # delete the temp files
    for fipath in new_fipath_list:
        os.remove(fipath)


if __name__ == '__main__':
    opt_parser = OptionParser()
    opt_parser.add_option('--cell', help='the folder including RNA-seq bed files of cell')
    opt_parser.add_option('--refgenes', help='the ref genes file path')
    options, args = opt_parser.parse_args()
    cell = options.cell
    refgene_file = options.refgenes
    output_read_count(cell, refgene_file)