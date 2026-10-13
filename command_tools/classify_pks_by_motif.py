#!/usr/bin/env python3
# coding=utf-8
"""Classify peaks by presence of a given motif."""
import argparse

from import_MAmotif import classify_MAnorm_peaks_by_motif


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='peak file path')
    parser.add_argument('-M', dest='motifscan', required=True, help='motifscan result file path')
    parser.add_argument('-n', dest='motifname', required=True, help='motif name')
    options = parser.parse_args(argv)
    classify_MAnorm_peaks_by_motif(options.pk, options.motifscan, options.motifname)


if __name__ == '__main__':
    command()
