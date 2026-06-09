"""Tests for tl.integer_cn — integer copy number + purity/ploidy from segments."""
import numpy as np
import pytest
import spatial_omics as sc


def _purity05_profile(a, seed=0):
    """3N gain -> obs 2.5, 1N loss -> obs 1.5 at tumour purity 0.5."""
    rng = np.random.default_rng(seed)
    prof = np.full(a.n_vars, 2.0) + rng.normal(0, 0.02, a.n_vars)
    i8 = np.where(a.var["chr"].values == "chr8")[0]
    prof[i8[20:60]] = 2.5 + rng.normal(0, 0.02, 40)
    prof[i8[90:110]] = 1.5 + rng.normal(0, 0.02, 20)
    return prof


def _seg(a, profile):
    sc.tl.segment(a, profile=profile, penalty_scale=2.0)
    return a


def test_integer_cn_fits_purity_and_states():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    _seg(a, _purity05_profile(a))
    tbl = sc.tl.integer_cn(a)                 # auto-fit purity
    r = a.uns["spatial_omics_integer_cn"]
    assert r["resolvable"] is True
    assert r["purity"] == pytest.approx(0.5, abs=0.05)
    c8 = tbl[tbl.chr == "chr8"]
    assert int(c8[c8.start == 20_000_000]["integer_cn"].iloc[0]) == 3   # gain -> CN3
    assert int(c8[c8.start == 90_000_000]["integer_cn"].iloc[0]) == 1   # loss -> CN1
    assert "integer_cn" in a.var


def test_integer_cn_given_purity():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    _seg(a, _purity05_profile(a))
    tbl = sc.tl.integer_cn(a, purity=0.5)
    c8 = tbl[tbl.chr == "chr8"]
    assert int(c8[c8.start == 20_000_000]["integer_cn"].iloc[0]) == 3
    assert int(c8[c8.start == 90_000_000]["integer_cn"].iloc[0]) == 1
    assert a.uns["spatial_omics_integer_cn"]["purity_source"] == "given"


def test_integer_cn_flat_genome_guard():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    rng = np.random.default_rng(1)
    _seg(a, np.full(a.n_vars, 2.0) + rng.normal(0, 0.02, a.n_vars))
    sc.tl.integer_cn(a)
    r = a.uns["spatial_omics_integer_cn"]
    assert r["resolvable"] is False
    assert "unidentifiable" in r["verdict"]
    assert set(np.unique(a.var["integer_cn"].values)) <= {-1, 2}   # all diploid (or unsegmented)


def test_integer_cn_requires_segment():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    with pytest.raises(KeyError):
        sc.tl.integer_cn(a)
