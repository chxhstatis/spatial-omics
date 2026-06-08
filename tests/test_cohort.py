"""Test tl.cohort_compare — cross-sample bulk CNA with sex(chrX) confound control."""
import numpy as np
import spatial_omics as sc


def _prep(seed, male=False):
    a = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=seed)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    if male:                                   # mimic a male: halve chrX copy number
        isx = (a.var["chr"] == "chrX").values
        a.layers["copy_number"][:, isx] *= 0.5
    return a


def test_cohort_compare_sex_and_correlation():
    res = sc.tl.cohort_compare({"F1": _prep(1), "M1": _prep(2, male=True), "F2": _prep(3)})
    # sex inferred from chrX
    assert res["inferred_sex"] == {"F1": "female", "M1": "male", "F2": "female"}
    assert res["chrX_mean_log2"]["M1"] < -0.5
    # correlation matrix is symmetric, unit diagonal, and chrX-excluded (autosomal)
    C = res["correlation"]
    assert list(C.index) == ["F1", "M1", "F2"]
    assert np.allclose(np.diag(C.values), 1.0)
    assert np.allclose(C.values, C.values.T)
    # recurrence + loci tables present
    assert {"chr", "start", "n_gain", "n_loss", "mean_log2"} <= set(res["recurrence"].columns)
    assert "locus" in res["recurrent_loci"].columns


def test_cohort_compare_accepts_list():
    a1, a2 = _prep(1), _prep(2)
    a1.uns["sample"] = "A"; a2.uns["sample"] = "B"
    res = sc.tl.cohort_compare([a1, a2])
    assert res["samples"] == ["A", "B"]
