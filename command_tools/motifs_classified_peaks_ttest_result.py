#!/usr/bin/env python3
# coding=utf-8
"""Run motif-classified M-value t-test / rank-sum tests for all motifs."""
import argparse

from import_MAmotif import motif_classified_pks_ttest


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='peak file path')
    parser.add_argument('-M', dest='motif', required=True, help='motifscan result file path')
    parser.add_argument('-n', dest='negative', action='store_true', default=False,
                        help='Negate M-values before testing')
    parser.add_argument('-c', dest='correction', default='benjamin',
                        help='pvalue correction: benjamin (default) or bonferroni')
    options = parser.parse_args(argv)
    motif_classified_pks_ttest(
        options.pk, options.motif,
        negative=options.negative, correction_type=options.correction)


if __name__ == '__main__':
    command()
