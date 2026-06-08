# spatial-omics — handoff runbook (read me first)

You are picking up a Python package, **`spatial-omics`** (import name `spatial_omics`),
for analysing **sparse spatial DNA-seq** (DBiT 50×50 / slide-DNA-seq): a spot × genomic-bin
fragment matrix → de-novo clones, per-spot copy number, spatially-resolved significance,
and multi-modal (RNA+DNA) variance decomposition.

Unlike inferCNV / CopyKAT / Numbat (which infer CNV from RNA *expression*), this works on
**measured DNA reads**. It is AnnData-native and follows the scanpy `pp`/`tl`/`pl` convention.

## 0. What's here

```
src/spatial_omics/   io · datasets · pp · tl · pl · _ref(bundled hg38 bias tracks)
tests/               pytest suite (clone-recovery ARI guard, vardecomp guard)
examples/            quickstart.py (synthetic)  run_real_data.py (REAL data)  fetch_public_data.sh
docs/                mkdocs-material site (4 tutorials + auto API)
pyproject.toml · CI (.github/workflows/ci.yml)
```

## 1. Install & verify (do this first)

```bash
cd spatial-omics
pip install -e ".[dev]"        # numpy pandas scipy scikit-learn anndata matplotlib
pytest -q                      # expect: 7 passed
python examples/quickstart.py  # runs full pipeline on synthetic data -> examples/_figures/
```

If `pytest` passes and quickstart writes 4 PNGs, the port is intact. Synthetic clone
recovery is ARI ≈ 0.96 — the test suite asserts this, so a regression will fail loudly.

## 2. Run on the REAL data (this is why it's on the server)

The data lives here on the server as **Stage-2 outputs** (per sample):
`top50x50_fragments_with_cb.tsv.gz` + `matched_spot_barcodes.tsv`.
(These come from the upstream `spatial_omics_pipeline/` stages 1–2: FASTQ → fragments
→ DBiT spatial mapping. If only FASTQ exists, run stages 1–2 first.)

```bash
python examples/run_real_data.py \
    --fragments /path/to/520_520/top50x50_fragments_with_cb.tsv.gz \
    --matched   /path/to/520_520/matched_spot_barcodes.tsv \
    --sample    520_520 \
    --outdir    results_520_520
```

This writes `results_520_520/520_520.h5ad` + clone / copy-number / significance figures.
Bias tracks (GC / mappability / blacklist, hg38 1 Mb) are **bundled** in `src/spatial_omics/_ref/`,
so no download is needed. Or in Python:

```python
import spatial_omics as so
adata = so.io.from_pipeline("…/top50x50_fragments_with_cb.tsv.gz", "…/matched_spot_barcodes.tsv", sample="520_520")
so.pp.load_reference_tracks(adata)        # bundled hg38 tracks
so.pp.bin_qc(adata); so.pp.correct_bias(adata); so.pp.normalize(adata)
so.tl.dual_smooth(adata); so.tl.call_clones(adata); so.tl.copy_number(adata)
so.tl.permutation_significance(adata, regions=["chr8","chr17","chr18"])
so.pl.spatial_clones(adata)
```

### Multi-modal (if matched spatial RNA exists)
Provide an RNA `AnnData` on the same grid (`obs['x_id']`, `obs['y_id']`, `obsm['spatial']`,
and a cell-type annotation `obs['cell_type']`, e.g. from RCTD), then:
```python
so.tl.variance_decomposition(adata_dna, adata_rna, tumor_types=("Tumor",), immune_types=("Immune",))
so.pl.variance_decomposition(adata_rna)
```

## 3. Known limits (do NOT over-claim — same caveats as the source pipeline)

- ~0.4× coverage → **CNA / clones only, no SNV** (KRAS point mutation is invisible).
- Low tumour purity → muted CNA amplitudes; bulk/clone CNA is the robust deliverable,
  per-spot subclones are weaker.
- `normalize` spatially densifies before building the pseudo-normal reference — required,
  because at this sparsity the raw cross-spot median is 0 for almost every bin.

## 4. Provenance

The three core methods (PC+xy double smoothing, permutation significance, de-novo clones)
and the bias-correction order are ported faithfully from the project's `spatial_omics_pipeline/`
stage3 code, themselves following slide-DNA-seq (Zhao et al. *Nature* 2022). Full ecosystem
plan and competitive analysis: see the project's `docs/ECOSYSTEM_PLAN.md`. Placeholders
(`your-org`, author, institution) in pyproject/LICENSE/CITATION should be filled before publishing.
