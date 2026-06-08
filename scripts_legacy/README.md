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
| `stage5_register.py` | `tl.register` |
| `ref/` data, `run_stage3.sh`, `README_stage3.md` | bundled `src/spatial_omics/_ref/` + docs site |

## Still here — utility + roadmap (the port backlog is now empty)

| File | Status |
|---|---|
| `stage3_cna_clone/compute_bin_gc.py` | **Utility** — builds the bundled GC reference from a genome FASTA (needed when adding a new genome, e.g. mm10 / T2T). Not a per-sample algorithm, so it stays a script. |
| `stage4_rna_integration/PLAN.md` | **Roadmap** — `tl.variance_decomposition` (needs matched RNA). |
| `stage5_wsi_registration/PLAN.md` | **Roadmap** — registration (done: `tl.register`) → annotation → anchoring (next). |
| `stage6_control_compare/PLAN.md` | **Roadmap** — control / Panel-of-Normals comparison. |

Every stage3–6 *algorithm* is now in the API (one implementation). What remains is a
reference-data builder + roadmap PLANs for features that need new data (RNA, controls).
