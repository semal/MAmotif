# coding=utf-8
# 将一组peaks按照peak是否与基因的promoter区域有overlap(重叠)将其分成两类
from optparse import OptionParser
import os
from import_MAmotif import classify_MAnorm_pks_by_promoter


def command():
    # -------- command -----------------------
    optparser = OptionParser()
    optparser.add_option('-p', dest='pk', help='pk file path')
    optparser.add_option('-r', dest='refgene', help='refgene file path')
    options, args = optparser.parse_args()
    pk = options.pk
    if pk.endswith(os.sep):
        pk = pk[:-1]
    refgene = options.refgene
    classify_MAnorm_pks_by_promoter(pk, refgene)


def test():
    from constant import mm9_refgenes_file
    os.chdir('/mnt/MAmotif/2.Processing/4.Esb4_and_Ese14_H3K4me3_MAmotif/5.MotifPeaks')
    for pk in os.listdir('.'):
        if not pk.endswith('.xls'):
            continue
        if 'Ese14' not in pk:
            continue
        cmd = \
            'python ~/MAmotif/command_tools/classify_pks_by_promoter.py -p %s -r %s &' % \
            (pk, mm9_refgenes_file)
        # classify_MAnorm_pks_by_promoter(pk, mm9_refgenes_file)
        print cmd


if __name__ == '__main__':
    # test()
    command()