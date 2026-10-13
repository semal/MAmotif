"""Smoke tests for the Python 3 MAmotif port."""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from MAmotif_pkg import (
    MAnormPeak,
    MAnormPeakSet,
    MAnormPeaksClassifier,
    Feature,
    classify_MAnorm_pks_by_promoter,
    read_MAnorm_peaks,
    read_refgenes,
    correct_pvalues,
)
from MAmotif import build_parser, main
import constant


DATA = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
PEAK_FILE = os.path.join(DATA, 'H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls')
HG19_REF = os.path.join(DATA, 'modified_hg19_refgene.txt')


def test_import_package():
    import MAmotif_pkg  # noqa: F401
    assert hasattr(MAmotif_pkg, 'read_MAnorm_peaks')


def test_read_manorm_peaks():
    peaks = read_MAnorm_peaks(PEAK_FILE)
    assert len(peaks) > 100
    assert peaks[0].mvalue is not None
    assert peaks[0].chrm.lower().startswith('chr')


def test_read_refgenes_sample():
    # Full hg19 file is large; only check that the reader returns Gene objects.
    genes = read_refgenes(HG19_REF)
    assert len(genes) > 1000
    g = genes[0]
    assert g.name
    assert g.chrm.startswith('chr')


def test_classify_by_promoter(tmp_path):
    # Work in temp dir because classify writes sibling output files.
    import shutil
    peak_copy = tmp_path / 'sample_peak_MAvalues.xls'
    shutil.copy(PEAK_FILE, peak_copy)
    # Use a tiny subset of peaks for speed: rewrite first 200 data lines
    with open(PEAK_FILE, encoding='utf-8') as src:
        lines = src.readlines()
    header = [ln for ln in lines if ln.startswith('#')]
    data = [ln for ln in lines if not ln.startswith('#')][:200]
    peak_copy.write_text(''.join(header + data), encoding='utf-8')

    # Also need a small refgene set overlapping those chroms — use full file once;
    # classification on 200 peaks is acceptable.
    promoter, distal = classify_MAnorm_pks_by_promoter(str(peak_copy), HG19_REF)
    assert os.path.exists(promoter)
    assert os.path.exists(distal)
    assert os.path.getsize(promoter) > 0


def test_ttest_on_synthetic_peaks():
    yes = MAnormPeakSet([
        _peak_with_m('chr1', 100, 200, 2.0),
        _peak_with_m('chr1', 300, 400, 2.5),
        _peak_with_m('chr1', 500, 600, 3.0),
    ])
    no = MAnormPeakSet([
        _peak_with_m('chr1', 700, 800, 0.1),
        _peak_with_m('chr1', 900, 1000, -0.2),
        _peak_with_m('chr1', 1100, 1200, 0.0),
    ])
    clf = MAnormPeaksClassifier(MAnormPeakSet(), Feature())
    clf.set_feature(yes, no)
    t_stat, t_p = clf.ttest_feature_classified_peaks()
    r_stat, r_p = clf.ranksum_feature_classified_peaks()
    assert t_stat is not None and t_p is not None
    assert t_p < 0.05
    assert r_stat is not None and r_p is not None


def test_correct_pvalues_benjamin():
    pvals = [0.001, 0.01, 0.2, None]
    corrected = correct_pvalues(pvals, 'benjamin')
    assert corrected[0] < corrected[1]
    assert corrected[3] is None
    assert all(c is None or c <= 1 for c in corrected)


def test_cli_help():
    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(['-h'])
    assert exc.value.code == 0


def test_cli_requires_args():
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2


def test_constant_paths():
    assert os.path.isdir(constant.data_folder)
    assert os.path.isfile(constant.hg19_refgenes_file)
    motifs = constant.jaspar3_motifs
    assert len(motifs) == 196


def _peak_with_m(chrm, start, end, mvalue):
    pk = MAnormPeak(chrm, start, end, summit=(end - start) // 2)
    pk.set_mvalue(mvalue)
    return pk
