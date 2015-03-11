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


if __name__ == '__main__':
    command()