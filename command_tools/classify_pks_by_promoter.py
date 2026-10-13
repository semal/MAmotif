#!/usr/bin/env python3
# coding=utf-8
"""Classify peaks by overlap with gene promoter regions."""
import argparse
import os

from import_MAmotif import classify_MAnorm_pks_by_promoter


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='peak file path')
    parser.add_argument('-r', dest='refgene', required=True, help='refgene file path')
    options = parser.parse_args(argv)
    pk = options.pk.rstrip(os.sep)
    classify_MAnorm_pks_by_promoter(pk, options.refgene)


if __name__ == '__main__':
    command()
