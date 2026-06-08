"""Tests for the rigor layer — confound controls must (a) confirm real structure on
synthetic data with planted clones, and (b) expose its own machinery honestly."""
import numpy as np
import spatial_omics as sc


def _prep(n_clones=2, coverage=0.4, seed=0):
    a = sc.datasets.simulate(coverage=coverage, n_clones=n_clones, seed=seed)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    return a


def test_morans_i_detects_planted_structure():
    a = _prep()
    I, p = sc.tl.morans_i(a.obs["cnv_burden"].values, a.obs["x_id"].values,
                          a.obs["y_id"].values, n_perm=199)
    assert np.isfinite(I) and I > 0.1          # planted clones are spatially autocorrelated
    assert p < 0.05


def test_spatial_heterogeneity_reports_verdict():
    a = _prep()
    res = sc.tl.spatial_heterogeneity(a, n_perm=199)
    assert {"moran_coverage", "moran_cnv_burden", "cna_exceeds_coverage"} <= set(res)
    assert isinstance(res["cna_exceeds_coverage"], bool)
    assert "spatial_omics_heterogeneity" in a.uns


def test_clone_diagnostics_runs_and_flags():
    a = _prep()
    res = sc.tl.clone_diagnostics(a)
    assert 0.0 <= res["clone_cna_distinctness"] <= 1.0
    assert isinstance(res["clones_likely_artifact"], bool)
    assert "spatial_omics_clone_diagnostics" in a.uns


def test_detect_channel_stripes_no_false_positive_on_random():
    a = _prep()
    sc.tl.permutation_significance(a, regions=["chr8"])
    res = sc.tl.detect_channel_stripes(a, "z_chr8")
    assert {"row_dominance", "col_dominance", "stripe_artifact"} <= set(res)
    # synthetic data has no microfluidic channel banding
    assert res["stripe_artifact"] is False
