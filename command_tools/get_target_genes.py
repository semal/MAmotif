#!/usr/bin/env python3
# coding=utf-8
"""Output target genes for a peak file given a RefGene annotation."""
import argparse

from import_MAmotif import output_peaks_target_genes


def command(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pk', help='peak file path')
    parser.add_argument('refgene', help='refgene annotation file')
    options = parser.parse_args(argv)
    output_peaks_target_genes(options.pk, options.refgene)


if __name__ == '__main__':
    command()
