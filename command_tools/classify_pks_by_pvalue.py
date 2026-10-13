#!/usr/bin/env python3
# coding=utf-8
"""Classify peaks by MAnorm p-value significance."""
import argparse
import os

from import_MAmotif import classify_MAnorm_pks_by_pvalue


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='peak file path')
    options = parser.parse_args(argv)
    classify_MAnorm_pks_by_pvalue(options.pk.rstrip(os.sep))


if __name__ == '__main__':
    command()
