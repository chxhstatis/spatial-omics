# Reproducing slide-DNA-seq (Nature 2022)

The strongest trust signal for a new tool is reproducing a landmark result. The
slide-DNA-seq paper (Zhao, Chiang et al., *Nature* 601:85–91, 2022) released raw data
under SRA **`PRJNA768453`** (mouse metastasis models, human brain tumour). Its official
analysis code is MATLAB, unmaintained since 2021, and not installable — here the same
result comes from a few lines of pip-installable Python.

!!! info "This tutorial needs a data download + the upstream pipeline"
    Unlike the synthetic tutorials, this one consumes real reads. The download is large;
    run the upstream steps on a machine with disk + bandwidth (e.g. your HPC). Once you
    have a `spot × bin` matrix, the downstream analysis below is identical to every other
    tutorial.

## 1. Download

```bash
# edit examples/fetch_public_data.sh to list the slide-DNA-seq SRR runs for PRJNA768453
# (SRA Run Selector: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA768453)
bash examples/fetch_public_data.sh public_data/slide_dna_seq
```

## 2. Reads → spot × bin matrix (upstream `workflow/`)

The bead-array assay differs from a DBiT chip, so the barcode-decoding step uses the
paper's bead barcodes rather than `workflow/stage2_spatial_map`. The output is the same
standard object: a spot × 1 Mb-bin fragment matrix with per-spot coordinates. Load it:

```python
import spatial_omics as sc

adata = sc.io.from_matrix(
    counts,                 # spots × bins fragment counts (np.ndarray / sparse)
    var=bins,               # DataFrame: chr, start, end
    obs=spots,              # DataFrame: x_id, y_id, total_frags (+ spatial coords)
    sample="mouse_met_1",
)
```

## 3. Downstream — identical to every other tutorial

```python
sc.pp.load_reference_tracks(adata)        # bundled hg38 GC / mappability / blacklist
sc.pp.bin_qc(adata); sc.pp.correct_bias(adata); sc.pp.normalize(adata)
sc.tl.dual_smooth(adata)                  # paper's pc_scores_smo_both
sc.tl.call_clones(adata)                  # paper Fig 2b/3b de-novo clones
sc.tl.copy_number(adata)
sc.tl.permutation_significance(adata, regions=["chr8", "chr17", "chr18"])  # paper Fig 2c
```

## 4. Confirm with the rigor layer, then compare to the figures

Because these samples are higher-purity than desmoplastic PDAC, the guards should
**confirm** the clones (unlike the low-purity negatives elsewhere):

```python
sc.tl.clone_diagnostics(adata)            # expect clones_likely_artifact = False
sc.pl.spatial_clones(adata)               # compare to paper Fig 2b / 3b
sc.pl.significance(adata, region="chr8")  # compare to paper Fig 2c
```

Recovering the paper's spatially-restricted clones and region significance — from an
installable Python package, validated by the same artifact guards used everywhere else —
is the correctness demonstration this tool is built to make.

## References

- Paper: <https://www.nature.com/articles/s41586-021-04217-4>
- Data: SRA `PRJNA768453`
- Original (MATLAB) code: <https://github.com/buenrostrolab/slide_dna_seq_analysis>
