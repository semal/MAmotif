# coding=utf-8
import os
from import_MAmotif import other_pks_classified_pks_ttest


def test():
    # ------------ test ------------------------
    import os
    os.chdir('F:\\MAmotif.py\\0. MAmotif_used_data')
    pk_all = 'H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues.xls'
    pk_promoter = 'H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls'
    pk_distal = 'H1hesc_H3K27ac_Broad_Rep2_stringent_distal_peak_MAvalues.xls'
    # classify_MAnorm_pks_by_other_peaks(pk_promoter, 'key_peaks\\human_all_tfs')
    # classify_MAnorm_pks_by_other_peaks(pk_distal, 'key_peaks\\human_all_tfs')
    other_pks_classified_pks_ttest(pk_all, 'key_peaks\\human_all_tfs')


def command():
    # -------- command -----------------------
    from optparse import OptionParser
    optparser = OptionParser()
    optparser.add_option('-p', dest='pk', help='pk file path')
    optparser.add_option('-b', dest='bind', help='other chipseq peak file')
    options, args = optparser.parse_args()
    pk = options.pk
    if pk.endswith(os.sep):
        pk = pk[:-1]
    bind_pks = options.bind
    other_pks_classified_pks_ttest(pk, bind_pks)


if __name__ == '__main__':
    command()
    # test()