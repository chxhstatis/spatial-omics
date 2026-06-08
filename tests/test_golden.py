"""Golden-output regression tests — pin known outputs so unintended numerical drift fails CI.

Tiers:
  * EXACT for numpy-deterministic quantities (simulate counts, planted clones).
  * TOLERANCE for image processing (he_purity) and modal-anchored copy number.
  * BEHAVIOURAL for the rigor guards (the differentiator): real clones confirmed,
    clone-free data rejected, channel stripes flagged. If one of these flips on a
    library bump, you WANT to know — investigate the guard's calibration.

Baselines captured 2026-06-08 (numpy-deterministic `simulate`).
"""
import numpy as np
import pytest
import spatial_omics as sc
from sklearn.metrics import adjusted_rand_score


@pytest.fixture(scope="module")
def pipe2():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    return a


def test_golden_simulate_deterministic():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    assert a.shape == (2500, 3044)
    assert int(np.asarray(a.layers["counts"]).sum()) == 3_025_044
    assert dict(a.obs["true_clone"].value_counts()) == {"Normal": 1688, "Clone1": 406, "Clone2": 406}


def test_golden_pipeline_metrics(pipe2):
    assert adjusted_rand_score(pipe2.obs["true_clone"], pipe2.obs["clone"]) > 0.90
    cn = np.asarray(pipe2.layers["copy_number"])
    cnf = cn[np.isfinite(cn)]
    assert np.median(cnf) == pytest.approx(2.0, abs=0.02)   # modal anchoring -> 2N
    assert np.mean(cnf) == pytest.approx(2.0, abs=0.05)


def test_golden_rigor_behaviour(pipe2):
    het = sc.tl.spatial_heterogeneity(pipe2, n_perm=99)
    assert het["moran_coverage"] < 0.10           # random tissue coverage ~ 0
    assert het["moran_cnv_burden"] > 0.50         # planted clones autocorrelate
    cd = sc.tl.clone_diagnostics(pipe2)
    assert cd["clones_likely_artifact"] is False  # real clones -> confirmed
    pipe2.obs["stripe"] = np.sin(pipe2.obs["y_id"].values * 1.3) * 3
    st = sc.tl.detect_channel_stripes(pipe2, "stripe")
    assert st["col_dominance"] > 0.90 and st["stripe_artifact"] is True


def test_golden_guard_rejects_false_clones():
    """The headline guarantee: with NO planted clones, naive clustering still reports
    clones, but clone_diagnostics rejects them."""
    a = sc.datasets.simulate(coverage=0.4, n_clones=1, seed=0)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    cd = sc.tl.clone_diagnostics(a)
    assert cd["clone_cna_distinctness"] < 0.05
    assert cd["clones_likely_artifact"] is True


def test_golden_he_purity_numbers():
    cv2 = pytest.importorskip("cv2")
    pytest.importorskip("skimage")
    import os
    import tempfile
    img = np.full((300, 400, 3), 245, np.uint8)
    img[60:240, 80:320] = (200, 150, 200)
    img[90:210, 120:280] = (150, 60, 150)
    cv2.rectangle(img, (80, 60), (320, 240), (40, 40, 220), 4)
    p = os.path.join(tempfile.mkdtemp(), "golden_he.jpg")
    cv2.imwrite(p, img)
    hp = sc.tl.he_purity(p)
    assert hp["box_color"] == "red"
    assert hp["nuclei_fraction"] == pytest.approx(0.424, abs=0.02)
    assert hp["stroma_fraction"] == pytest.approx(0.071, abs=0.02)
