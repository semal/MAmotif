# coding=utf-8
# MAnorm的结果给每个M值都对应了一个pvalue用来衡量变化的显著性，p值越小说明变化越显著，利用给定的p值将一组peaks分成两类
from optparse import OptionParser
import os
from import_MAmotif import classify_MAnorm_pks_by_pvalue


def command():
    optparser = OptionParser()
    optparser.add_option('-p', dest='pk', help='pk file path')
    options, args = optparser.parse_args()
    pk = options.pk
    if pk.endswith(os.sep):
        pk = pk[:-1]
    classify_MAnorm_pks_by_pvalue(pk)

if __name__ == '__main__':
    command()