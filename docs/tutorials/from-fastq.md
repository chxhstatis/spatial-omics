# From FASTQ (upstream)

`spatial_omics` analyses a **spot × bin matrix**. Producing that matrix from raw reads
is the upstream pipeline (the "Cell Ranger" layer): linker/ME trimming, spatial
barcode extraction, alignment, dedup, spatial mapping to the DBiT grid.

That pipeline lives in [`spatial_omics_pipeline/`](https://github.com/your-org/spatial_omics)
(stages 1–2) and is containerised for HPC/cloud. Outline:

```text
FASTQ
  └─ stage1  linker filter → strip ME → extract BC2+BC1 → trim R1 → bowtie2(hg38)
             → dedup → q30 → per-barcode fragments
  └─ stage2  DBiT 50×50: split BC2/BC1 (16 bp) → match 192 refs (≤1 mm)
             → crop 50×50 → real coordinates
  → top50x50_fragments_with_cb.tsv.gz  +  matched_spot_barcodes.tsv
```

Then hand the two Stage-2 files to spatial_omics:

```python
import spatial_omics as sc
adata = sc.io.from_pipeline(
    "top50x50_fragments_with_cb.tsv.gz",
    "matched_spot_barcodes.tsv",
    sample="520_520",
)
# attach bias tracks (GC / mappability / blacklist) — bundled tracks used by default
sc.pp.load_reference_tracks(adata)                    # packaged hg38 tracks
# sc.pp.load_reference_tracks(adata, "path/to/custom/ref")   # or your own
```

!!! warning "Coverage limits"
    At ~0.4× genome coverage the robust deliverables are **CNA and clones**.
    Point mutations (e.g. KRAS) are invisible — that needs targeted/deep sequencing.

*This page is a stub — full upstream container instructions land in P5.*
