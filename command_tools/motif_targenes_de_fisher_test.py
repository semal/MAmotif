#!/usr/bin/env python3
# coding=utf-8
"""Fisher exact test of motif target genes vs differentially expressed genes."""
import argparse

from import_MAmotif import motif_targenes_de_fisher_test


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', dest='pk', required=True, help='MAnorm peak file path')
    parser.add_argument('-M', dest='ms', required=True, help='MotifScan result file path')
    parser.add_argument('-r', dest='ref', required=True, help='RefSeq gene file path')
    parser.add_argument('-d', dest='de', required=True, help='Differentially expressed gene list')
    options = parser.parse_args(argv)
    motif_targenes_de_fisher_test(options.pk, options.ms, options.de, options.ref)


if __name__ == '__main__':
    command()
