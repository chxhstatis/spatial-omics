# scripts_legacy — what's here (and what isn't)

After the duplication audit (2026-06-08), the stage3–6 algorithms that the package now
implements were **removed from here** so that each algorithm has exactly ONE
implementation (the `src/spatial_omics` API). The originals are preserved in this repo's
git history and in the project archive (`99_归档Archive/旧版代码/spatial_omics_pipeline/`).

## Now in the package API (removed from scripts_legacy)

| Old script | Canonical API |
|---|---|
| `stage3_per_spot_cnv.py` | `tl.copy_number` |
| `stage3_spatial_clones.py` | `tl.dual_smooth` + `tl.call_clones` + `tl.permutation_significance` |
| `stage3_spatial_heterogeneity.py` | `tl.spatial_heterogeneity` / `tl.morans_i` |
| `pick_normal_spots.py` | `pp.pick_normal_spots` |
| `stage3_bulk_cna.py` | `pp.normal_anchor` + `tl.copy_number` (bulk = pseudobulk) |
| `compare_samples_cna.py` | `tl.cohort_compare` |
| `stage5_purity_he.py` | `tl.he_purity` |
| `ref/` data, `run_stage3.sh`, `README_stage3.md` | bundled `src/spatial_omics/_ref/` + docs site |

## Still here — NOT yet in the API (the port backlog) + roadmap

| File | Status |
|---|---|
| `stage3_cna_clone/compute_bin_gc.py` | **Utility** — builds the bundled GC reference from a genome FASTA (needed when adding a new genome, e.g. mm10 / T2T). Keep. |
| `stage5_wsi_registration/stage5_register.py` (+ `README_stage5.md`) | **Port target** — H&E ROI-box → 50×50 grid registration. Becomes `tl.register` (uses the `[he]` extra, like `tl.he_purity`). |
| `stage4_rna_integration/PLAN.md` | **Roadmap** — `tl.variance_decomposition` (needs matched RNA). |
| `stage5_wsi_registration/PLAN.md` | **Roadmap** — registration → annotation → anchoring. |
| `stage6_control_compare/PLAN.md` | **Roadmap** — control / Panel-of-Normals comparison. |

When a port-target lands in the API, delete it here too (keep "one implementation").
