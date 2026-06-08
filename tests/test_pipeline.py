"""End-to-end tests on synthetic data with known ground truth.

These also serve as the regression guard for the ported stage3 algorithms: if a
change breaks bias correction, smoothing or clustering, the ARI / amplitude
assertions fail.
"""
import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

import spatial_omics as sc


@pytest.fixture(scope="module")
def processed():
    adata = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    sc.pp.bin_qc(adata)
    sc.pp.correct_bias(adata)
    sc.pp.normalize(adata)
    sc.tl.dual_smooth(adata)
    sc.tl.call_clones(adata, kmin=2, kmax=6)
    sc.tl.copy_number(adata)
    return adata


def test_simulate_shape_and_truth():
    adata = sc.datasets.simulate(coverage=0.5, n_clones=2, seed=1)
    assert adata.n_obs > 0 and adata.n_vars > 2000
    assert {"gc", "mappability", "blacklist"} <= set(adata.var.columns)
    assert set(adata.obs["true_clone"]) == {"Normal", "Clone1", "Clone2"}
    assert "spatial" in adata.obsm


def test_relative_is_finite_after_normalize():
    adata = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=2)
    sc.pp.bin_qc(adata)
    sc.pp.correct_bias(adata)
    sc.pp.normalize(adata)
    keep = adata.var["pass_qc"].values
    rel = np.asarray(adata.layers["relative"])[:, keep]
    # densification must keep the reference well-defined (the sparse-median bug)
    assert np.isfinite(rel).mean() > 0.95


def test_clone_recovery(processed):
    ari = adjusted_rand_score(processed.obs["true_clone"], processed.obs["clone"])
    assert ari > 0.7, f"clone recovery degraded (ARI={ari:.3f})"
    assert "Normal" in set(processed.obs["clone"])


def test_copy_number_amplitude(processed):
    cn = np.asarray(processed.layers["copy_number"])
    chr8 = (processed.var["chr"] == "chr8").values
    tum = (processed.obs["clone"] != "Normal").values
    nrm = (processed.obs["clone"] == "Normal").values
    # chr8 carries a simulated gain -> tumour CN should exceed normal
    assert np.nanmean(cn[tum][:, chr8]) > np.nanmean(cn[nrm][:, chr8])


def test_permutation_significance(processed):
    sc.tl.permutation_significance(processed, regions=["chr8"], n_perm=30, seed=0)
    assert "z_chr8" in processed.obs
    assert np.isfinite(processed.obs["z_chr8"]).any()


def test_no_depth_confound(processed):
    r = processed.uns["spatial_omics_dual_smooth"]["depth_confound_r_pc1"]
    assert abs(r) < 0.6, f"PC1 confounded by coverage (r={r:.2f})"


def test_variance_decomposition_recovers_drivers(processed):
    rna = sc.datasets.simulate_rna(processed, n_genes=150,
                                   frac_genetic=0.3, frac_microenv=0.3, seed=1)
    sc.tl.variance_decomposition(processed, rna)
    v = rna.var
    gen = v[v.true_driver == "genetic"]
    mic = v[v.true_driver == "microenv"]
    assert gen["vd_genetic"].mean() > gen["vd_microenv"].mean()
    assert mic["vd_microenv"].mean() > mic["vd_genetic"].mean()
    assert (gen["vd_dominant"] == "genetic").mean() > 0.8
    assert (mic["vd_dominant"] == "microenv").mean() > 0.8
