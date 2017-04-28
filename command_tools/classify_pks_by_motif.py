# coding=utf-8
# 将一组peaks(peak是定义的一个组蛋白修饰区域，本质上是一个DNA序列)按照是否存在某个moitf分成两类
from optparse import OptionParser
from import_MAmotif import classify_MAnorm_peaks_by_motif


def command():
    # --------------------command--------------------------------------------
    optparser = OptionParser()
    optparser.add_option('-p', dest='pk', help='peak file path')
    optparser.add_option('-M', dest='motifscan', help='motifscan result file path')
    optparser.add_option('-n', dest='motifname', help='motif name')
    options, args = optparser.parse_args()
    pk = options.pk
    motifscan = options.motifscan
    motifname = options.motifname
    classify_MAnorm_peaks_by_motif(pk, motifscan, motifname)
    # ----------------------------------------------------------------


if __name__ == '__main__':
    command()
    pass