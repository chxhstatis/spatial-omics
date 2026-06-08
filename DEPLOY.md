# Deploying the documentation site

The docs are a **static** site built by `mkdocs-material` (the source is `docs/` +
`mkdocs.yml`). A static site needs no server — host it on a free CDN. Two routes; pick
by audience.

## Local preview

```bash
pip install -e ".[docs]"
mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build            # static site -> ./site
```

## Before first deploy: replace the placeholders

```bash
# set your GitHub org/repo and (optional) custom domain
grep -rl "your-org" mkdocs.yml pyproject.toml README.md | xargs sed -i '' 's/your-org/YOUR_GH_ORG/g'
```

(Also set `authors` email in `pyproject.toml` / `CITATION.cff` if you want it public.)

## Route A — GitHub Pages (default; best for an international audience)

Already wired in CI (`.github/workflows/ci.yml`, the `docs` job): every push to `main`
runs `mkdocs gh-deploy --force`, which builds and publishes to the `gh-pages` branch.

1. Push the repo to GitHub.
2. Settings → Pages → Source = `gh-pages` branch.
3. Site goes live at `https://YOUR_GH_ORG.github.io/spatial-omics/`.
4. Custom domain: add a `CNAME` (Settings → Pages) and a DNS `CNAME` record.

Free, zero-maintenance, global CDN. **Caveat: GitHub Pages can be slow or intermittently
blocked from mainland China.**

## Route B — Cloudflare Pages (recommended if your users are in China)

Cloudflare Pages is free, has a China-accessible CDN, and builds the same static site.
No code change needed — configure in the dashboard:

1. Push the repo to GitHub (or GitLab).
2. Cloudflare dashboard → **Workers & Pages → Create → Pages → Connect to Git** → select the repo.
3. Build settings:
   - **Framework preset:** None
   - **Build command:** `pip install -e ".[docs]" && mkdocs build`
   - **Build output directory:** `site`
   - **Environment variable:** `PYTHON_VERSION = 3.11`
4. Deploy. Every push to `main` rebuilds automatically.
5. Custom domain: Pages → Custom domains → add yours (Cloudflare manages the DNS + TLS).

> For best mainland-China latency, use a `.com`/`.org` domain proxied through Cloudflare;
> a domain with mainland ICP filing (备案) on a domestic CDN is the most robust but needs
> filing. Cloudflare Pages on a global domain is usually good enough and needs no filing.

You can run both routes at once (GitHub Pages international + Cloudflare Pages for China)
from the same repo — they don't conflict.

## What you do NOT need

- **A server / Mac mini for the docs** — it is static; a CDN serves it. Keep the Mac mini
  for CI, an internal mirror, or lab compute.
- **A backend** — there is no online analysis here. Users `pip install spatial-omics` and
  run locally. (An optional upload-and-analyze web app would be a separate, later project
  on on-demand cloud compute — not part of this docs site.)
