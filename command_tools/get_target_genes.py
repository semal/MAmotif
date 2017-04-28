# coding=utf-8
# 这个脚本用来输出给定peaks的靶基因
import os
from import_MAmotif import *


def get_target_genes(pk_folder, refgene_fp):
    for fn in os.listdir(pk_folder):
        fp = os.sep.join([pk_folder, fn])
        output_peaks_target_genes(fp, refgene_fp)


def test():
    os.chdir('/mnt/MAmotif/3.Analysis/2.mESCImportantGenesRelated/Esb4_H3K4me3_target_genes')
    refgene_fp = '/home/zhaohui/MAmotif/data/modified_mm9_refgene.txt'
    pk_fd = '/mnt/MAmotif/2.Processing/4.Esb4_and_Ese14_H3K4me3_MAmotif/5.MotifPeaks/promoter'
    # get_target_genes(pk_fd, refgene_fp)
    for pk_name in os.listdir(pk_fd):
        if not pk_name.endswith('.xls'):
            continue
        print 'python ~/MAmotif/command_tools/get_target_genes.py %s %s &' % \
              (os.sep.join([pk_fd, pk_name]), refgene_fp)


def command():
    # ----------command-----------------------------
    if '-h' in sys.argv:
        print 'python get_target_genes.py pk refgene'
        exit(1)
    pk1 = sys.argv[1]
    ref = sys.argv[2]
    output_peaks_target_genes(pk1, ref)


if __name__ == '__main__':
    command()
    # test()