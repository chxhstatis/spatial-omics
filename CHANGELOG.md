## [0.2.0] - 2026-06-08
### Consolidated
- Merged the standalone `spatial_omics_pipeline` (stage scripts) into this one canonical repo:
  stage1-2 -> `workflow/`, stage3-6 -> `scripts_legacy/` (reference implementations).
- Authors set; single bundled reference data.
### Added — rigor layer (`tl.rigor`)
- `spatial_heterogeneity`, `clone_diagnostics`, `detect_channel_stripes`, `morans_i`:
  confound controls / artifact guards that reject coverage/smoothing/channel-stripe
  false clones. 4 new tests; full suite green.

# Changelog

## [0.1.0] — unreleased

First scaffold. AnnData-native package porting the `spatial_omics_pipeline` stage3
algorithms into a scanpy-style API.

### Added
- `io.from_pipeline` / `io.from_matrix` — build a standard AnnData (spot × bin).
- `datasets.simulate` — synthetic spatial DNA-seq with ground-truth clones.
- `pp.bin_qc`, `pp.correct_bias`, `pp.normalize` (spatial densification + pseudo-normal).
- `tl.dual_smooth` (PC + xy double smoothing), `tl.call_clones` (CH-selected k),
  `tl.copy_number`, `tl.permutation_significance`.
- `pl` spatial maps + clone profiles.
- Test suite (clone-recovery ARI guard), runnable `examples/quickstart.py`, mkdocs
  site, GitHub Actions CI.
