#!/usr/bin/env python3
# coding=utf-8
"""MAmotif pipeline: integrate MAnorm M-values with MotifScan results."""
from __future__ import annotations

import argparse
import os
import sys

from MAmotif_pkg import classify_MAnorm_pks_by_promoter, motif_classified_pks_ttest


def run_MAmotif_pipeline(
        comparison_pk, motifscan_result,
        refgene_file='', correction_type='benjamin', neg=False, no_fd=True):
    comparison_pk = os.path.abspath(comparison_pk)
    motifscan_result = os.path.abspath(motifscan_result)
    MAnorm_name = comparison_pk.split(os.sep)[-2]
    if not no_fd:
        if not os.path.exists(MAnorm_name):
            os.mkdir(MAnorm_name)
        os.chdir(MAnorm_name)
    # classify pk into promoter_pk and distal-pk
    if refgene_file == '':
        motif_classified_pks_ttest(
            comparison_pk, motifscan_result, correction_type=correction_type, negative=neg)
    else:
        promoter_pk, distal_pk = classify_MAnorm_pks_by_promoter(comparison_pk, refgene_file)
        # motifs classified peaks test for pk, promoter_pk and distal_pk
        motif_classified_pks_ttest(
            comparison_pk, motifscan_result, correction_type=correction_type, negative=neg)
        motif_classified_pks_ttest(
            promoter_pk, motifscan_result, correction_type=correction_type, negative=neg)
        motif_classified_pks_ttest(
            distal_pk, motifscan_result, correction_type=correction_type, negative=neg)


def build_parser():
    parser = argparse.ArgumentParser(
        description='Integrate MotifScan results with MAnorm binding-signal changes.')
    parser.add_argument('-p', dest='pk', required=True,
                        help='MAnorm comparison peak file path')
    parser.add_argument('-r', dest='refgene', default='',
                        help='RefGene annotation file (optional; enables promoter/distal split)')
    parser.add_argument('-M', dest='motif', required=True,
                        help='MotifScan result file path')
    parser.add_argument('-n', dest='negative', action='store_true', default=False,
                        help='Negate M-values before testing (swap comparison direction)')
    parser.add_argument('-c', dest='correction', default='benjamin',
                        help='p-value correction: benjamin (default), bonferroni, or other=none')
    parser.add_argument('-d', dest='nfd', action='store_true', default=False,
                        help='Do not create a MAnorm comparison output folder')
    return parser


def command(argv=None):
    options = build_parser().parse_args(argv)
    run_MAmotif_pipeline(
        options.pk,
        options.motif,
        options.refgene,
        correction_type=options.correction,
        neg=options.negative,
        no_fd=options.nfd,
    )


def main(argv=None):
    command(argv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
