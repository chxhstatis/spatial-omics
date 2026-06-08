# Artifact guards: don't report false clones

On sparse, low-purity tissue, naive clustering finds "clones" **even when there are
none** — driven by tissue-coverage density, the smoothing kernel, or microfluidic
channel stripes. `spatial_omics` ships the confound controls that reject these before
you report them. This is the honest core of the toolkit, and the demo that matters
most. (Runnable script: `examples/artifact_guards.py`.)

```python
import numpy as np
import spatial_omics as sc

def run(n_clones, seed=1):
    a = sc.datasets.simulate(coverage=0.4, n_clones=n_clones, seed=seed)
    sc.pp.bin_qc(a); sc.pp.correct_bias(a); sc.pp.normalize(a)
    sc.tl.dual_smooth(a); sc.tl.call_clones(a); sc.tl.copy_number(a)
    return a
```

## 1. `clone_diagnostics` — are the clones real?

`call_clones` will always return *something*. The guard asks whether those clusters
differ at chromosome-arm scale, and whether `k` was chosen on merit (not pinned to the
search ceiling).

```python
for label, n in [("REAL clones", 2), ("NO clones", 1)]:
    a = run(n)
    cd = sc.tl.clone_diagnostics(a)
    print(label, "-> k =", cd["chosen_k"],
          "| distinctness =", round(cd["clone_cna_distinctness"], 3),
          "| artifact =", cd["clones_likely_artifact"])
```

```text
REAL clones -> k = 3 | distinctness = 0.084 | artifact = False
NO clones   -> k = 2 | distinctness = 0.036 | artifact = True
```

Even with **no planted clones**, naive `call_clones` reports two — but
`clone_diagnostics` flags them: their copy-number profiles are near-identical
(distinctness below the floor), so they are a coverage/smoothing artifact, not biology.

!!! danger "This is the failure mode that sinks papers"
    inferCNV / CopyKAT / the official slide-DNA-seq MATLAB code have no such guard. On
    low-purity samples they will happily output a clone map. `clone_diagnostics` is what
    keeps your customers' results honest and publishable.

## 2. `detect_channel_stripes` — is it a microfluidic artifact?

DBiT coverage is dominated by per-channel efficiency. A set of weak X- or Y-channels
shows up as full rows/columns in a per-spot map (a region z-score, coverage, …) and
masquerades as spatial biology.

```python
a = run(2)
y = a.obs["y_id"].values
a.obs["channel_striped"] = np.sin(y * 1.3) * 3        # a fake Y-channel banding
a.obs["random_map"]      = np.random.default_rng(1).normal(0, 1, a.n_obs)

print(sc.tl.detect_channel_stripes(a, "channel_striped"))  # stripe_artifact=True, Y-channels
print(sc.tl.detect_channel_stripes(a, "random_map"))       # stripe_artifact=False
```

```text
channel_striped -> stripe_artifact=True  (col_dominance=0.99, axis='Y-channels (columns)')
random_map      -> stripe_artifact=False (col_dominance=0.02)
```

This is exactly how the deep-sequencing report caught a "significant chr17 deletion"
that was really a column of low-efficiency channels.

## 3. `spatial_heterogeneity` — is structure above the coverage baseline?

Tissue density is itself spatially autocorrelated. Real spatial CNA heterogeneity must
exceed it. We compute Moran's I of coverage and of the CNV burden, with a true
random-permutation null.

```python
a = run(2)
het = sc.tl.spatial_heterogeneity(a, n_perm=999)
print(het["moran_coverage"], het["moran_cnv_burden"], het["cna_exceeds_coverage"])
```

!!! warning "Use the unsmoothed burden for the honest test"
    `cnv_burden` from `tl.copy_number` with `k_spatial>1` is **smoothing-inflated** —
    neighbouring spots share input, so their burden autocorrelates by construction. For
    the strict verdict, compute an unsmoothed burden first:

    ```python
    sc.tl.copy_number(a, k_spatial=1)      # no spatial densification
    het = sc.tl.spatial_heterogeneity(a)   # now the verdict is honest
    ```

    On real low-purity sections this is the step that exposes coverage/smoothing-driven
    structure as *not* clonal — the negative result that protects you from over-claiming.

## The point

The same three functions **confirm** planted clones on clean data and **reject** the
artifacts that appear on real low-purity samples. Run them before reporting any spatial
clone. Honest negatives ("no robust clones here, and here is why") are a feature, not a
failure — they are what makes the platform's results trustworthy.

→ Next: [Reproducing slide-DNA-seq (Nature 2022)](reproduce-nature2022.md).
