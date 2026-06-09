# Methods roadmap

What to add to grow `spatial-omics` from a CNA/clone core (v0.1.0) into a real
spatial-genomics method ecosystem — organized by **data dependency** (the binding
constraint) and **priority**.

## Principle: demand-driven, with a guard

Ecosystems like scanpy/Seurat grew one function at a time, each solving a real analysis
that was stuck — not by pre-building breadth that nobody uses (dead code is the main way
academic tools rot). So: **add a method when an analysis needs it.** And per this project's
core, **every new method ships with a guard / honest failure mode** — it must say "not
trustworthy here" when the data can't support it. Never add a method that quietly
over-claims.

## Where v0.1.0 already is

`io` · `datasets.simulate` · `pp`(bin_qc / correct_bias / normalize / normal_anchor /
pick_normal_spots) · `tl`(dual_smooth / call_clones / copy_number / permutation_significance
/ **rigor**: spatial_heterogeneity / clone_diagnostics / detect_channel_stripes / morans_i /
he_purity / cohort_compare / register) · `pl`(maps + clone profiles).

This is a **bin-level CNA + de-novo clone + rigor** core. The gaps below are what an
ecosystem needs next.

---

## Tier A — feasible NOW (no new data) · build first

| Method | What & why | API | Priority | Status |
|---|---|---|---|---|
| **CNA segmentation** | CBS / HMM / changepoint → smooth copy-number **segments** instead of bin-level values (cf. Ginkgo, HMMcopy, DNAcopy). The current `copy_number` is bins + median filter; segmentation is the foundation for integer CN, breakpoints, and focal-event calling. | `tl.segment` | **P0** | todo |
| **Integer copy number + ploidy** | Assign integer CN states to segments and estimate ploidy/scale (depth-ratio-based; no BAF at this depth). Makes CNA calls publication-grade and is the input clone phylogeny needs. | `tl.integer_cn` | **P0** (after segment) | todo |
| **Purity from sequencing** | Estimate tumour fraction from the CNA-amplitude distribution (the spread of relative ratios), independent of H&E. Makes the "low purity" argument quantitative and drives normal anchoring. Pairs with `he_purity`. | `tl.estimate_purity` | P1 | todo |
| **Spatial domains** | Spatial-neighbour graph + cluster CNA into spatially-coherent **domains** (Leiden on a spatial+CNA graph; or via squidpy). "Regions of shared CNA" — more honest than per-spot clones at low resolution. | `tl.spatial_domains` | P1 | todo |
| **Clone spatial geometry** | Boundary / interface analysis, spatial mixing index, clone adjacency. Spatial arrangement is the whole point of the assay — quantify it. | `tl.clone_geometry` | P1 | todo |
| **Recurrent-CNA scoring** | Cohort recurrent-region scoring with significance (GISTIC-like), extending `cohort_compare`; must exclude centromere/chrX artifacts (already learned). | `tl.recurrent_cna` | P2 | partial (cohort_compare) |
| **Benchmark vs RNA-CNV** | Run inferCNV/CopyKAT-style RNA-inferred CNV on matched data and report concordance — the trust/validation story (we infer from *measured DNA*, they from expression). | `examples/benchmark` | P2 | todo |
| **Plotting breadth** | Genome-wide heatmap (spots × bins), domain maps, recurrence track, (later) clone phylogeny — as `pl` functions. | `pl.*` | P1–P2 | incremental |

> **Recommended next build:** `tl.segment` → `tl.integer_cn`. Everything downstream
> (phylogeny, focal events, dosage) needs real segments, and it lifts CNA quality now.

## Tier B — needs same-section RNA (the multimodal soul) · highest value

| Method | What & why | API | Priority | Data |
|---|---|---|---|---|
| **RNA load + joint object** | Read matched spatial RNA on the same 50×50 grid; align to the DNA object (MuData / shared `.obsm`). | `io.from_rna` | P0-when-RNA | same-section RNA |
| **Cell-type deconvolution** | Thin wrapper over RCTD / cell2location → per-spot tumour / immune / stroma proportions = the microenvironment-density covariates. | `tl.deconvolve` | P0-when-RNA | RNA |
| **★ Variance decomposition** | `expression ~ subclone + tumour_density + immune_density` per gene; partition variance into **intrinsic (genetic) vs extrinsic (microenvironment)**. This is Zhao 2022 Fig 4 and the single method that separates this tool from every RNA-inferred-CNV tool. | `tl.variance_decomposition` | **P0-when-RNA (flagship)** | RNA |
| **RNA-defined normal anchor** | Use RNA-defined stromal/immune spots as the external normal for `pp.normal_anchor` — solves the low-purity "no internal normal" problem we proved. | feeds `pp.normal_anchor` | P0-when-RNA | RNA |
| **CNA dosage analysis** | Does CNA gain/loss drive expression up/down (cis-dosage)? Links genotype → phenotype. | `tl.cna_dosage` | P1-when-RNA | RNA |
| **RNA orientation anchor** | Use RNA tissue structure to resolve the H&E registration orientation (the unresolved 8-fold ambiguity in `tl.register`). | feeds `tl.register` | P1-when-RNA | RNA |

## Tier C — needs controls / external normal

| Method | What & why | API | Priority | Data |
|---|---|---|---|---|
| **Panel-of-Normals** | Build a PoN from control samples; anchor against it. The correct normal baseline for low-purity tumours (internal anchoring proven to fail). | `pp.build_pon` / `normal_anchor(pon=…)` | P0-when-controls | control samples |
| **Differential CNA** | Tumour-vs-control differential copy number with significance. | `tl.differential_cna` | P1-when-controls | controls |

## Tier D — needs higher purity / experimental change (deferred)

| Method | Why deferred |
|---|---|
| **Clone phylogeny / evolution** (MEDICC2-style from CNA) | Needs ≥2 distinct-CNA clones as tree nodes; current low-purity PDAC has none. Unlock via purity enrichment (microdissection / tumour-rich selection). |
| **SNV / point mutations** | Out of scope at ~0.4–4× depth; needs targeted deep sequencing (separate assay). |

## Tier E — interop / infrastructure (no data dep, broadens the ecosystem)

| Method | What & why | API | Priority |
|---|---|---|---|
| **More platform readers** | Visium HD, Stereo-seq, slide-DNA-seq bead arrays, 10x → one analysis core for many platforms. | `io.read_*` | P2 |
| **SpatialData / squidpy interop** | Make objects SpatialData-compatible; reuse squidpy spatial tools. | `io` | P2 |
| **Multi-genome** | mm10 / T2T support — generalise `compute_bin_gc`, bundle versioned references with checksums. | `pp` / `_ref` | P2 |
| **Upstream containerisation** | nextflow / snakemake + Docker for FASTQ→matrix — the platform's "Cell Ranger" (one-command, HPC/cloud). | `workflow/` | P1 (platform value) |
| **More bias covariates** | Replication-timing track, etc. (rep/Tn5 are GC-redundant — add only if it measurably helps). | `pp` | P3 |

---

## What to build, in order

1. **Now (Tier A):** `tl.segment` + `tl.integer_cn` — real CNA calling, the base for everything. Then `tl.spatial_domains` + `tl.clone_geometry` + plotting.
2. **Reserve the interfaces** for Tier B so RNA plugs in instantly: a `tl.variance_decomposition` stub + the joint-object I/O contract.
3. **On data arrival:** RNA → Tier B (flagship variance decomposition); controls → Tier C (PoN); purity enrichment → Tier D (phylogeny).
4. **In parallel (platform):** containerise `workflow/` (Tier E) — direct adoption value for chip customers.

Each method: one `pp`/`tl`/`pl` function, scanpy-style signature, a numpy docstring, a
test (functional + a golden baseline), and a guard. See `CONTRIBUTING.md`.
