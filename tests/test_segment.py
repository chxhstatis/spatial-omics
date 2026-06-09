"""Tests for tl.segment — piecewise-constant CNA segmentation."""
import numpy as np
import pytest
import spatial_omics as sc


def _planted_profile(a, seed=0):
    rng = np.random.default_rng(seed)
    prof = np.full(a.n_vars, 2.0) + rng.normal(0, 0.03, a.n_vars)
    i8 = np.where(a.var["chr"].values == "chr8")[0]
    prof[i8[20:60]] = 2.6 + rng.normal(0, 0.03, 40)   # gain block
    prof[i8[90:110]] = 1.5 + rng.normal(0, 0.03, 20)   # loss block
    return prof


def test_segment_recovers_planted_steps():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    tbl = sc.tl.segment(a, profile=_planted_profile(a), penalty_scale=2.0)
    res = a.uns["spatial_omics_segments"]
    assert res["flat_genome"] is False
    assert res["var_explained"] > 0.8
    c8 = tbl[tbl.chr == "chr8"]
    # golden: exact boundaries of the planted blocks, correct segment means
    gain = c8[(c8.start == 20_000_000) & (c8.end == 60_000_000)]
    loss = c8[(c8.start == 90_000_000) & (c8.end == 110_000_000)]
    assert len(gain) == 1 and gain["mean"].iloc[0] == pytest.approx(2.6, abs=0.05)
    assert len(loss) == 1 and loss["mean"].iloc[0] == pytest.approx(1.5, abs=0.05)


def test_segment_flat_genome_guard():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    rng = np.random.default_rng(1)
    flat = np.full(a.n_vars, 2.0) + rng.normal(0, 0.03, a.n_vars)   # no events
    sc.tl.segment(a, profile=flat, penalty_scale=2.0)
    res = a.uns["spatial_omics_segments"]
    assert res["flat_genome"] is True
    assert res["n_segments"] <= 23 + 1   # ~one per chromosome (autosomes + chrX)


def test_segment_writes_var_and_runs_on_pseudobulk():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.copy_number(a)
    tbl = sc.tl.segment(a)                 # pseudobulk of copy_number
    assert "segment" in a.var and "seg_mean" in a.var
    assert {"chr", "start", "end", "n_bins", "mean"} <= set(tbl.columns)
