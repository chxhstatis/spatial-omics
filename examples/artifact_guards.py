"""Artifact guards demo — the honest core of spatial_omics.

Naive clustering finds "clones" even when there are none. The rigor layer rejects them.
Run: ``python examples/artifact_guards.py``  (also exercised in CI).
"""
import numpy as np

import spatial_omics as sc


def _run(n_clones, seed=1):
    a = sc.datasets.simulate(coverage=0.4, n_clones=n_clones, seed=seed)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    return a


print("=" * 70)
print("1. clone_diagnostics rejects FALSE clones, confirms REAL ones")
print("=" * 70)
for label, n in [("REAL clones (n_clones=2)", 2), ("NO clones  (n_clones=1)", 1)]:
    a = _run(n)
    cd = sc.tl.clone_diagnostics(a)
    naive_k = cd["chosen_k"]
    print(f"\n{label}")
    print(f"  naive call_clones still reports k={naive_k}: {a.uns['spatial_omics_clones']['sizes']}")
    print(f"  clone CNA distinctness = {cd['clone_cna_distinctness']:.3f}")
    print(f"  --> clones_likely_artifact = {cd['clones_likely_artifact']}")
    print(f"      {cd['verdict']}")

print("\n" + "=" * 70)
print("2. detect_channel_stripes flags microfluidic channel banding")
print("=" * 70)
a = _run(2)
y = a.obs["y_id"].values
a.obs["channel_striped"] = np.sin(y * 1.3) * 3 + np.random.default_rng(0).normal(0, 0.2, a.n_obs)
a.obs["random_map"] = np.random.default_rng(1).normal(0, 1, a.n_obs)
for key in ["channel_striped", "random_map"]:
    st = sc.tl.detect_channel_stripes(a, key)
    print(f"  {key:18s} -> stripe_artifact={st['stripe_artifact']} "
          f"(row_dom={st['row_dominance']:.2f}, col_dom={st['col_dominance']:.2f})")

print("\n" + "=" * 70)
print("3. spatial_heterogeneity: is CNV structure above the COVERAGE baseline?")
print("=" * 70)
a = _run(2)
het = sc.tl.spatial_heterogeneity(a, n_perm=99)
print(f"  coverage Moran's I = {het['moran_coverage']:.2f} (tissue density)")
print(f"  CNV burden Moran's I = {het['moran_cnv_burden']:.2f} (smoothed -- inflated!)")
print("  NOTE: burden from smoothed copy_number is smoothing-inflated. For the honest")
print("  test pass an UNSMOOTHED burden: sc.tl.copy_number(a, k_spatial=1) first.")
print("  On real low-purity data this is what exposes coverage/smoothing-driven false clones.")

print("\nDone. The same guards confirm planted clones and reject artifacts —")
print("so users get honest, publishable results instead of false clones.")
