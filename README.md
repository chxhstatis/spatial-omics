# spatial_omics

[![CI](https://github.com/chxhstatis/spatial-omics/actions/workflows/ci.yml/badge.svg)](https://github.com/chxhstatis/spatial-omics/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Spatial copy-number & clone analysis for sparse spatial DNA-seq** (DBiT / slide-DNA-seq).

Unlike inferCNV / CopyKAT / Numbat — which infer copy number from RNA *expression* —
`spatial_omics` works on **measured DNA reads**: a spot × genomic-bin fragment matrix. It is
AnnData-native and follows the scanpy `pp`/`tl`/`pl` convention, so it slots into the
scverse ecosystem (scanpy, squidpy, SpatialData).

```python
import spatial_omics as sc

adata = sc.datasets.simulate(coverage=0.4, n_clones=2)   # or sc.io.from_pipeline(...)
sc.pp.bin_qc(adata); sc.pp.correct_bias(adata); sc.pp.normalize(adata)
sc.tl.dual_smooth(adata)        # PC + xy double smoothing — beats sparsity
sc.tl.call_clones(adata)        # de-novo clones (Calinski–Harabasz selects k)
sc.tl.copy_number(adata)        # per-spot copy number
sc.tl.permutation_significance(adata, regions=["chr8", "chr18"])
sc.pl.spatial_clones(adata)
```

## Install

```bash
pip install -e ".[dev]"     # from source (PyPI release planned)
python examples/quickstart.py
pytest
```

## What's inside

| Module | Purpose |
|---|---|
| `io` | `from_pipeline` / `from_matrix` → standard AnnData object |
| `datasets` | `simulate()` synthetic ground-truth data; public-data loaders (planned) |
| `pp` | bin QC, GC/mappability bias correction, normalize; **`normal_anchor`** (+ signal-collapse guard) and **`pick_normal_spots`** |
| `tl` | **dual_smooth**, **call_clones**, **copy_number**, **permutation_significance**; **`he_purity`** (H&E cellularity, `[he]` extra), **`cohort_compare`** (cross-sample), **`register`** (H&E↔grid, `[he]` extra) |
| `tl.rigor` | **confound controls / artifact guards** — `spatial_heterogeneity`, `clone_diagnostics`, `detect_channel_stripes`, `morans_i` |
| `pl` | spatial clone / copy-number / significance maps, clone CNA profiles |

The methods that make sparse spatial genomics work (double smoothing, layered bias
correction, permutation significance) are carried over from slide-DNA-seq (Zhao et al.
*Nature* 2022) and the `spatial_omics_pipeline` stage3 code.

### Why the rigor layer matters

On sparse, low-purity tissue, naive clustering manufactures "clones" out of tissue-coverage
density, the smoothing kernel, and microfluidic channel stripes. Before reporting a clone,
`tl.rigor` proves it is real: `spatial_heterogeneity` asks whether CNV autocorrelation
exceeds the *coverage* baseline (Moran's I, true permutation null); `clone_diagnostics` asks
whether clones differ at *arm scale* and whether `k` was chosen on merit rather than pinned
to the search ceiling; `detect_channel_stripes` flags grid-aligned channel banding. The same
functions confirm planted clones on synthetic data and reject the artifacts that appear on
real low-purity samples — users get honest, publishable results instead of false clones.

Repository: `src/spatial_omics/` is the analysis package; `workflow/` is the upstream
FASTQ→matrix step; `scripts_legacy/` holds the stage3–6 reference implementations. See
`MERGE_MAP.md` for how the previous codebases were consolidated.

## Status

Alpha (v0.1.0). Validated on synthetic data: clone recovery ARI ≈ 0.96. Roadmap and
the broader ecosystem plan: see `docs/` and the project `ECOSYSTEM_PLAN.md`.

## Citing

See [`CITATION.cff`](CITATION.cff). A methods preprint is planned.

## License

MIT.
