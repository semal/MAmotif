# coding=utf-8
# 这个脚本用来输出给定peaks的靶基因
from import_MAmotif import *


def get_target_genes(pk_folder, refgene_fp):
    for fn in os.listdir(pk_folder):
        fp = os.sep.join([pk_folder, fn])
        output_peaks_target_genes(fp, refgene_fp)


def test():
    os.chdir('F:\\MAmotif.py\\7. primed_naive_hESC\\primed_hESC_H3K4me3_genes')
    refgene_fp = 'F:\\MAmotif_src\\data\modified_hg19_refgene.txt'
    get_target_genes('peaks', refgene_fp)


def command():
    # ----------command-----------------------------
    if '-h' in sys.argv:
        print 'python get_target_genes.py pk refgene'
        exit(1)
    pk1 = sys.argv[1]
    ref = sys.argv[2]
    output_peaks_target_genes(pk1, ref)


if __name__ == '__main__':
    # command()
    test()