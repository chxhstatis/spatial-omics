# Changelog

All notable changes to `spatial-omics`. Format: [Keep a Changelog](https://keepachangelog.com).

## [Unreleased]

## [0.1.0] — 2026-06-08

First public release. An AnnData-native, scanpy-style toolkit for **sparse spatial
DNA-seq** (DBiT / slide-DNA-seq): a spot × genomic-bin fragment matrix → de-novo clones,
per-spot copy number, and spatially-resolved significance — with confound controls that
reject the false clones naive methods produce on low-purity tissue.

### Data & I/O
- `io.from_pipeline` / `io.from_matrix` / `read_h5ad` — build a standard AnnData (spot × bin).
- `datasets.simulate` — synthetic spatial DNA-seq with ground-truth clones (no patient data needed).
- Bundled hg38 reference tracks (GC / mappability / blacklist) — works out of the box.

### Preprocessing (`pp`)
- `bin_qc`, `correct_bias` (GC + mappability), `normalize` (spatial densification + pseudo-normal).
- `normal_anchor` (+ signal-collapse guard) and `pick_normal_spots` — external/data-driven
  normal anchoring with a guard that detects when no genuine internal normal exists.

### Tools (`tl`)
- `dual_smooth` (PC + xy double smoothing), `call_clones` (Calinski–Harabasz selects k),
  `copy_number`, `permutation_significance`.
- **Rigor layer** (`tl.rigor`): `spatial_heterogeneity` (CNV Moran vs coverage baseline,
  true-permutation null), `clone_diagnostics` (CNA distinctness + CH-boundary),
  `detect_channel_stripes` (microfluidic banding), `morans_i`.
- `he_purity` — independent H&E colour-deconvolution cellularity/stroma estimate (`[he]` extra).
- `cohort_compare` — cross-sample bulk CNA; sample similarity on autosomes only (chrX
  excluded so samples don't cluster by sex); inferred sex; recurrence + recurrent-loci tables.

### Plotting (`pl`)
- Spatial clone / copy-number / significance maps; clone CNA profiles.

### Project
- `workflow/` (upstream FASTQ→matrix) and `scripts_legacy/` (stage3–6 reference implementations).
- 16-test suite (clone-recovery ARI guard + rigor/anchor/he/cohort), runnable
  `examples/quickstart.py` and `examples/artifact_guards.py`, three mkdocs tutorials,
  GitHub Actions CI (test + docs deploy).
