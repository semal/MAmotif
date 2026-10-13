#!/usr/bin/env python3
# coding=utf-8
"""Classify peaks by overlap with another ChIP-seq peak set, then test M-values."""
import argparse
import os

from import_MAmotif import other_pks_classified_pks_ttest


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='peak file path')
    parser.add_argument('-b', dest='bind', required=True, help='other ChIP-seq peak folder/file')
    options = parser.parse_args(argv)
    other_pks_classified_pks_ttest(options.pk.rstrip(os.sep), options.bind)


if __name__ == '__main__':
    command()
