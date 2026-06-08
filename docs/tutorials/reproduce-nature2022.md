# Reproducing slide-DNA-seq (Nature 2022)

The strongest trust signal for a new tool is reproducing a landmark result. The
slide-DNA-seq paper (Zhao, Chiang et al., *Nature* 2021/2022, "Spatial genomics
enables multi-modal study of clonal heterogeneity in tissues") released raw data
under **SRA `PRJNA768453`** — mouse metastasis models and human brain tumour.

The plan for this flagship tutorial (P3 on the roadmap):

1. Download a slide-DNA-seq array from `PRJNA768453`.
2. Run the upstream pipeline → spot × bin matrix → `sc.io.from_pipeline`.
3. `pp.bin_qc → correct_bias → normalize`, then `tl.dual_smooth → call_clones`.
4. Compare the recovered clone map and per-region significance to the paper's
   Figs 2–3, computed here with **a few lines of pip-installable Python** instead
   of the original MATLAB.

!!! info "Why this matters"
    The official analysis code (`buenrostrolab/slide_dna_seq_analysis`) is MATLAB,
    unmaintained since 2021, and not installable. Reproducing its results in
    spatial_omics demonstrates correctness *and* that the field finally has an open,
    scriptable tool for this data type.

*This page is a stub — it is filled in once the public data is processed (P3).*

## References

- Paper: <https://www.nature.com/articles/s41586-021-04217-4>
- Data: SRA `PRJNA768453`
- Original code: <https://github.com/buenrostrolab/slide_dna_seq_analysis>
