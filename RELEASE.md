# Release checklist — publishing `spatial-omics` to PyPI

The name **`spatial-omics`** is free on PyPI (verified 2026-06-08). Import name is
`spatial_omics`. Below is the exact path from this repo to `pip install spatial-omics`.

## 0. Pre-flight (already verified locally)

- ✅ `python -m build` produces `dist/spatial_omics-0.1.0.tar.gz` + `…-py3-none-any.whl`.
- ✅ Bundled `_ref/` data (GC / mappability / blacklist) is included in the wheel.
- ✅ `twine check dist/*` → PASSED (both).
- ✅ Clean-install test: the installed wheel imports, loads its bundled reference, and runs
  end-to-end (`pip install --no-deps --target … dist/*.whl`).
- ✅ 16 tests pass.

## 1. Fix the placeholders (one-time, before first upload)

```bash
# set your GitHub org everywhere
grep -rl "your-org" pyproject.toml mkdocs.yml README.md | xargs sed -i '' 's/your-org/YOUR_GH_ORG/g'
# (optional) set a contact email in pyproject [project].authors and CITATION.cff
```

`[project.urls]` Homepage/Documentation/Issues should point at the real repo before upload.

## 2. Get API tokens

- PyPI account → Account settings → **API tokens** → create a token (scope: entire account
  for the first upload, then project-scoped). Same for **TestPyPI** (test.pypi.org).
- Store as `~/.pypirc` or pass via `TWINE_USERNAME=__token__ TWINE_PASSWORD=pypi-…`.

## 3. Build clean

```bash
rm -rf dist build *.egg-info
python -m build
twine check dist/*
```

## 4. Dry-run on TestPyPI first

```bash
twine upload --repository testpypi dist/*
# verify it installs from TestPyPI in a fresh venv:
python -m venv /tmp/v && /tmp/v/bin/pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ spatial-omics
/tmp/v/bin/python -c "import spatial_omics as sc; sc.datasets.simulate()"
```

## 5. Publish to PyPI

```bash
twine upload dist/*
# now anyone can:  pip install spatial-omics
```

## 6. Tag + GitHub release (→ Zenodo DOI)

```bash
git tag -a v0.1.0 -m "spatial-omics 0.1.0"
git push origin main --tags
```

- On GitHub, enable the **Zenodo–GitHub integration** (zenodo.org → linked accounts →
  toggle the repo on) **before** publishing the release. Then create a GitHub Release from
  the `v0.1.0` tag — Zenodo mints a DOI automatically. `.zenodo.json` (in this repo) supplies
  the metadata. Put the DOI badge in the README.

## 7. (Optional) Bioconda

Bioinformatics users expect conda. After the PyPI release, open a PR to
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) adding
`recipes/spatial-omics/meta.yaml` (a draft is in `conda-recipe/meta.yaml` — update the
`sha256` to the PyPI sdist hash, printed by `openssl dgst -sha256 dist/*.tar.gz`).

## Version bumps for later releases

Edit `version` in **both** `pyproject.toml` and `src/spatial_omics/__init__.py`
(`__version__`), add a `CHANGELOG.md` entry, then repeat steps 3–6.

## Notes

- Wheel is pure-Python (`py3-none-any`) — one wheel works on all platforms. No compilation.
- `[he]` extra (opencv + scikit-image) is optional and lazily imported, so the base install
  stays lightweight; `he_purity` raises a clear message if the extra is missing.
