# Contributing to spatial-omics

Thanks for helping maintain `spatial-omics`. This guide gets a new contributor productive
fast and keeps the codebase healthy as it grows.

## Philosophy (read this first)

1. **Honesty is the product.** On sparse, low-purity tissue, naive clustering invents
   clones. Our differentiator is the **rigor layer** (`tl.rigor`) that *rejects* those
   artifacts. Any new analysis method should ship with — or respect — a guard that says
   "this result is not trustworthy here" when it isn't. Never make a method that quietly
   over-claims.
2. **Stand on scverse.** AnnData-native, scanpy `pp` / `tl` / `pl` conventions. Don't invent
   a new object model; a new contributor who knows scanpy should feel at home immediately.
3. **One algorithm, one implementation.** If you port something from `scripts_legacy/`,
   delete the script (the original stays in git history). No parallel copies that can drift.

## Dev setup

```bash
git clone <repo> && cd spatial-omics
python -m pip install -e ".[dev,he]"   # dev = pytest/ruff/build; he = opencv+scikit-image
pytest -q                              # full suite (~5 min; uses synthetic data, no patient data)
ruff check src                         # lint
mkdocs serve                           # preview the docs at :8000
```

## The AnnData data-model contract

Every function reads/writes a shared object. Keep to these keys so functions compose:

| Slot | Key | Written by | Meaning |
|---|---|---|---|
| `.layers` | `counts` | `io` / `datasets` | spot × bin fragment counts (raw) |
| `.layers` | `corrected` | `pp.correct_bias` | GC + mappability corrected |
| `.layers` | `relative` | `pp.normalize` / `pp.normal_anchor` | divided by pseudo-normal |
| `.layers` | `log2_ratio` | `tl.dual_smooth` | smoothed log2 |
| `.layers` | `copy_number` | `tl.copy_number` | per-spot copy number (~2N) |
| `.var` | `chr,start,end,gc,mappability,blacklist,pass_qc` | `pp` | per-bin covariates + QC flag |
| `.obs` | `x_id,y_id,total_frags` | `io` / `datasets` | grid coords + coverage |
| `.obs` | `clone,cluster,cnv_burden,tissue_density` | `tl` | per-spot results |
| `.obsm` | `spatial,X_pca,X_pca_smoothed,X_he` | `tl` | coords / embeddings |
| `.uns` | `spatial_omics_*` | every tool | result summaries / diagnostics |

If you need a new intermediate, add a clearly-named layer/obs/uns key and document it here.

## Adding a method

1. Put it in the right module: `pp/` (preprocess), `tl/` (analysis), `pl/` (plots). One
   public function, scanpy-style signature `f(adata, *, ...)`, mutating `adata` and
   returning it (or a small result dict stored in `.uns`).
2. **numpy-style docstring** (mkdocstrings renders it into the API docs). State what keys it
   reads and writes.
3. Export it: add to the module's import block / `__all__`.
4. **Add a test** (see below). A method without a test will not be merged.
5. If it's an optional-dependency method (image/etc.), import the dep **lazily inside the
   function** and add the dep to a `[project.optional-dependencies]` extra in `pyproject.toml`
   (pattern: `tl.he_purity`, `tl.register` → `[he]`).
6. Add a one-line `CHANGELOG.md` entry under `[Unreleased]`.

## Testing

Two kinds, both required to stay green:

- **Functional** (`tests/test_*.py`): does it run, produce the right keys, basic correctness
  on synthetic data with known ground truth (e.g. clone-recovery ARI > 0.9).
- **Golden / regression** (`tests/test_golden.py`): pins known numeric outputs so
  *unintended* drift fails CI. Tiers: exact for numpy-deterministic quantities, tolerance
  for image/anchored values, **behavioural for the rigor guards** (real clones confirmed,
  clone-free data rejected).

**When you intentionally change an algorithm** and a golden value shifts: re-capture the
baseline (run the snippet at the top of `test_golden.py`'s capture history), update the
expected number, and say so in the PR + CHANGELOG. Never loosen a tolerance just to make a
red golden test pass — that defeats the point.

Tests use only `datasets.simulate` (synthetic) — never commit patient data or real images.

## Style

- `ruff` (line length 100, py39 target). Run `ruff check src` before pushing.
- Type hints on public signatures; `from __future__ import annotations` at the top.
- Keep functions focused; helpers prefixed `_`.

## Docs

`docs/` is mkdocs-material. Tutorials are runnable — the code in `examples/*.py` is executed
in CI, so tutorial code must stay in sync with it. API reference is auto-generated from
docstrings. See `DEPLOY.md` to publish.

## Release

See `RELEASE.md` (build → TestPyPI → PyPI → git tag → GitHub release/Zenodo → bioconda).
Bump `version` in **both** `pyproject.toml` and `src/spatial_omics/__init__.py`.

## Where the unported pieces live

`scripts_legacy/` holds only a reference-data builder (`compute_bin_gc.py`) and roadmap
PLANs (stage4 RNA variance-decomposition, stage5 annotation/anchoring, stage6 controls).
The originals of everything already in the API are in git history and the project archive.
`MERGE_MAP.md` records the consolidation; `scripts_legacy/README.md` maps old script → API.
