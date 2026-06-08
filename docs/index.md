# spatial_omics

**Spatial copy-number & clone analysis for sparse spatial DNA-seq** (DBiT / slide-DNA-seq).

Tools like inferCNV, CopyKAT and Numbat infer copy number from RNA *expression*.
`spatial_omics` is built for platforms that **measure DNA directly** — a spot × genomic-bin
fragment matrix — and turns it into de-novo clones, per-spot copy number, and
spatially-resolved significance maps.

It is **AnnData-native** and follows the scanpy `pp` / `tl` / `pl` convention, so it
drops straight into the scverse ecosystem (scanpy, squidpy, SpatialData).

```python
import spatial_omics as sc

adata = sc.datasets.simulate(coverage=0.4, n_clones=2)   # or sc.io.from_pipeline(...)
sc.pp.bin_qc(adata); sc.pp.correct_bias(adata); sc.pp.normalize(adata)
sc.tl.dual_smooth(adata)          # PC + xy double smoothing  (beats sparsity)
sc.tl.call_clones(adata)          # de-novo clones, CH selects k
sc.tl.copy_number(adata)          # per-spot integer-scaled copy number
sc.tl.permutation_significance(adata, regions=["chr8", "chr18"])
sc.pl.spatial_clones(adata)
```

## Why it exists

Sparse spatial genomics has three hard problems; spatial_omics packages a solution to each:

| Problem | Method | API |
|---|---|---|
| Each spot is too sparse to cluster | **Double smoothing** in PC-space *and* xy-space | `tl.dual_smooth` |
| Coverage ≠ copy number (GC / mappability / depth bias) | Layered quantile correction + pseudo-normal reference | `pp.correct_bias`, `pp.normalize` |
| "Is this amplification significant *here*?" | **Permutation** empirical-null z-scores per region | `tl.permutation_significance` |

## Status

Alpha. Validated on synthetic data (clone recovery ARI ≈ 0.96). See
[Getting started](tutorials/getting-started.md).
