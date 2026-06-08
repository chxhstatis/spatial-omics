# Getting started

This tutorial runs the whole pipeline on a synthetic sample with a **known**
clone layout, so you can see exactly what each step recovers. No data download
needed. (The runnable script is `examples/quickstart.py`.)

## 1. Get a dataset

```python
import spatial_omics as sc

adata = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
adata
# AnnData: n_obs (spots) × n_vars (1 Mb bins)
#   obs: x_id, y_id, total_frags, true_clone
#   var: chr, start, end, gc, mappability, blacklist
#   obsm: 'spatial'    layers: 'counts'
```

`simulate` models the real difficulties: ~0.4× genome coverage, GC/mappability
bias, and spatial tumour domains. On your own data use instead:

```python
adata = sc.io.from_pipeline("top50x50_fragments_with_cb.tsv.gz",
                            "matched_spot_barcodes.tsv", sample="520_520")
```

## 2. Preprocess

```python
sc.pp.bin_qc(adata)        # flag usable bins (GC, mappability, blacklist, coverage)
sc.pp.correct_bias(adata)  # GC + mappability quantile correction -> layers['corrected']
sc.pp.normalize(adata)     # spatial densify + library norm + pseudo-normal ref -> layers['relative']
```

!!! note "Why densify?"
    At ~0.4× coverage almost every (spot, bin) is 0, so the cross-spot median —
    the pseudo-normal reference — would be 0 everywhere. `normalize` first averages
    each spot over its nearest spatial neighbours (`spatial_k`), which is what makes
    sparse spatial DNA-seq analysable at all.

## 3. Analyse

```python
sc.tl.dual_smooth(adata)     # PC-space + xy-space double smoothing, then PCA
sc.tl.call_clones(adata)     # de-novo KMeans; Calinski–Harabasz picks k; flattest = Normal
sc.tl.copy_number(adata)     # per-spot copy number (anchored to 2N) -> layers['copy_number']
sc.tl.permutation_significance(adata, regions=["chr8", "chr18"])  # spatial z-score maps

adata.uns["spatial_omics_clones"]["sizes"]   # {'Normal': 1709, 'Clone1': 397, 'Clone2': 394}
```

## 4. Plot

```python
sc.pl.spatial_clones(adata)                      # clones painted on the grid
sc.pl.spatial_copy_number(adata, chrom="chr8")   # copy number map
sc.pl.significance(adata, region="chr8")         # signed permutation z-score
sc.pl.clone_profiles(adata)                      # genome-wide CNA per clone
```

## 5. Check it worked

Because the ground truth is known here, you can score the recovery:

```python
from sklearn.metrics import adjusted_rand_score
adjusted_rand_score(adata.obs["true_clone"], adata.obs["clone"])   # ≈ 0.96
```

Save the object for later (standard `.h5ad`, readable by scanpy/squidpy):

```python
adata.write("sim.h5ad")
```
