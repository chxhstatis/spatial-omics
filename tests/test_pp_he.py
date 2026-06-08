"""Tests for pp.normal_anchor / pp.pick_normal_spots and tl.he_purity."""
import numpy as np
import pytest
import spatial_omics as sc


def _prep():
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    return a


def test_pick_normal_spots_returns_mask():
    a = _prep()
    mask = sc.pp.pick_normal_spots(a, quantile=0.3)
    assert mask.dtype == bool and mask.shape == (a.n_obs,)
    assert 20 <= mask.sum() < a.n_obs


def test_normal_anchor_sets_guard():
    a = _prep()
    mask = sc.pp.pick_normal_spots(a, quantile=0.3)
    sc.pp.normal_anchor(a, mask)
    r = a.uns["spatial_omics_normal_anchor"]
    assert {"bulk_sd_before", "bulk_sd_after", "signal_collapsed", "n_reference_spots"} <= set(r)
    assert isinstance(r["signal_collapsed"], bool)
    assert r["n_reference_spots"] == int(mask.sum())


def test_he_purity_on_synthetic_image(tmp_path):
    cv2 = pytest.importorskip("cv2")
    pytest.importorskip("skimage")
    # synthetic H&E: white slide, red ROI rectangle, purple/pink "tissue" inside
    img = np.full((300, 400, 3), 245, np.uint8)            # white background (BGR)
    img[60:240, 80:320] = (200, 150, 200)                  # pink-ish tissue (BGR)
    img[90:210, 120:280] = (150, 60, 150)                  # denser purple nuclei region
    cv2.rectangle(img, (80, 60), (320, 240), (40, 40, 220), 4)   # red box (BGR)
    p = tmp_path / "synthetic_he.jpg"
    cv2.imwrite(str(p), img)
    res = sc.tl.he_purity(str(p))
    assert res["box_color"] == "red"
    assert 0.0 <= res["nuclei_fraction"] <= 1.0
    assert 0.0 <= res["stroma_fraction"] <= 1.0
    assert "H_over_E" in res


def test_register_on_synthetic_image(tmp_path):
    cv2 = pytest.importorskip("cv2")
    img = np.full((300, 400, 3), 245, np.uint8)
    img[60:240, 80:320] = (200, 150, 200)
    cv2.rectangle(img, (80, 60), (320, 240), (40, 40, 220), 4)   # red ROI box
    p = tmp_path / "he.jpg"; cv2.imwrite(str(p), img)
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
    res = sc.tl.register(a, str(p))
    assert res["box_color"] == "red"
    assert len(res["box_corners"]) == 4
    assert a.obsm["X_he"].shape == (a.n_obs, 2)          # per-spot pixel coords stored
    assert "tissue_density" in a.obs
    assert isinstance(res["orientation_resolved"], bool)
    assert len(res["all_orientation_r"]) == 8            # all 8 orientations scored
