"""Rigor layer — confound controls and artifact guards for sparse spatial DNA-seq.

These are the checks that separate *real* spatial copy-number structure from the
artifacts that naive clustering produces on sparse, low-purity data. They are the
honest core of the toolkit: before reporting a spatial clone, prove it is not just
tissue-coverage density, the smoothing kernel, or a microfluidic channel stripe.

Functions
---------
spatial_heterogeneity(adata)
    Moran's I of CNV burden vs the coverage (tissue-density) baseline, with a true
    random-permutation null. Verdict ``cna_exceeds_coverage`` is the honest test:
    real spatial CNA heterogeneity must exceed the coverage autocorrelation.
clone_diagnostics(adata)
    Are the de-novo clones real? Two guards: (1) clone CNA *distinctness* (do clusters
    differ at arm scale, or are they near-identical?), (2) CH model-selection boundary
    (was k chosen on merit, or pinned to the search ceiling kmax?).
detect_channel_stripes(adata, key)
    Flags grid-aligned (full row/column) banding in a per-spot map — the signature of
    microfluidic channel-efficiency artifacts that masquerade as spatial CNA.

Background: low tumour purity caps CNA amplitude; per-spot smoothing and tissue
coverage are both spatially autocorrelated and will fake "clones" if not controlled.
See the project reports and ``scripts_legacy/stage3_cna_clone`` for the derivations.
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------- Moran's I
def _edges8(x_id, y_id):
    """8-neighbour adjacency among present spots -> array of (i, j) index pairs (i<j)."""
    pos = {(int(x), int(y)): i for i, (x, y) in enumerate(zip(x_id, y_id))}
    edges = []
    for i, (x, y) in enumerate(zip(x_id, y_id)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                j = pos.get((int(x) + dx, int(y) + dy))
                if j is not None and i < j:
                    edges.append((i, j))
    return np.asarray(edges, dtype=int)


def morans_i(values, x_id, y_id, *, n_perm=999, seed=0):
    """Moran's I spatial autocorrelation with an 8-neighbour weight and a TRUE random
    label-permutation null. Returns ``(I, p_perm)``; ``(nan, nan)`` if degenerate.

    The permutation reassigns the observed values to spot positions at random, which
    destroys spatial arrangement while preserving the value distribution.
    """
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    if ok.sum() < 20:
        return np.nan, np.nan
    x_id = np.asarray(x_id)[ok]; y_id = np.asarray(y_id)[ok]; vv = v[ok]
    edges = _edges8(x_id, y_id)
    if len(edges) < 10 or np.std(vv) == 0:
        return np.nan, np.nan
    z = vv - vv.mean(); denom = (z ** 2).sum()
    I = (len(vv) / len(edges)) * ((z[edges[:, 0]] * z[edges[:, 1]]).sum() / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for p in range(n_perm):
        zp = rng.permutation(vv); zp = zp - zp.mean()
        null[p] = (len(vv) / len(edges)) * ((zp[edges[:, 0]] * zp[edges[:, 1]]).sum() / denom)
    p_val = (1 + int((null >= I).sum())) / (1 + n_perm)
    return float(I), float(p_val)


def spatial_heterogeneity(adata, *, burden_key="cnv_burden", coverage_key="total_frags",
                          n_perm=999, seed=0):
    """Is there spatial CNA heterogeneity beyond the tissue-coverage baseline?

    Computes Moran's I of (a) coverage = tissue density and (b) the per-spot CNV burden,
    and reports the honest verdict ``cna_exceeds_coverage``. Real spatial CNA structure
    must autocorrelate MORE than coverage does; if it does not, the apparent structure
    is tissue density (and, if ``cnv_burden`` came from a spatially-smoothed
    ``copy_number``, the smoothing kernel) rather than biology.

    .. note::
       ``cnv_burden`` from ``tl.copy_number`` with ``k_spatial>1`` is smoothing-inflated.
       For the strict test, also run ``tl.copy_number(adata, k_spatial=1)`` and pass its
       burden as ``burden_key`` — the unsmoothed value is the fair comparison to coverage.

    Stores the result in ``adata.uns['spatial_omics_heterogeneity']`` and returns it.
    """
    obs = adata.obs
    if burden_key not in obs:
        raise KeyError(f"{burden_key!r} not in adata.obs — run tl.copy_number first")
    mi_cov, p_cov = morans_i(obs[coverage_key].values, obs["x_id"].values, obs["y_id"].values,
                             n_perm=n_perm, seed=seed)
    mi_bur, p_bur = morans_i(obs[burden_key].values, obs["x_id"].values, obs["y_id"].values,
                             n_perm=n_perm, seed=seed)
    exceeds = bool(np.isfinite(mi_bur) and np.isfinite(mi_cov) and mi_bur > mi_cov)
    res = {"moran_coverage": mi_cov, "p_coverage": p_cov,
           "moran_cnv_burden": mi_bur, "p_cnv_burden": p_bur,
           "cna_exceeds_coverage": exceeds,
           "verdict": ("spatial CNA heterogeneity above coverage baseline" if exceeds
                       else "NO heterogeneity beyond tissue coverage / smoothing")}
    adata.uns["spatial_omics_heterogeneity"] = res
    return res


# ----------------------------------------------------------------------------- clone guards
def clone_diagnostics(adata, *, distinct_thresh=0.2, distinct_floor=0.05):
    """Are the de-novo clones real, or a clustering artifact?

    Two independent guards (run ``tl.call_clones`` first):

    1. **CNA distinctness** — fraction of bins where the across-clone spread of the
       median-centred clone profiles exceeds ``distinct_thresh``. Real subclones differ
       at chromosome-arm scale; an artifact has near-identical, flat profiles
       (distinctness ~ 0).
    2. **CH boundary** — if the Calinski-Harabasz-selected ``k`` equals the search
       ceiling ``kmax`` with no interior peak, ``k`` is a boundary artifact (CH has no
       over-splitting penalty), not evidence of that many clones.

    Sets ``adata.uns['spatial_omics_clone_diagnostics']`` with a
    ``clones_likely_artifact`` flag and returns it.
    """
    cl = adata.uns.get("spatial_omics_clones")
    if cl is None:
        raise KeyError("run tl.call_clones first")
    prof = cl["profiles"]  # bins x clones (log2-ish, median-centred)
    P = np.asarray(prof.values, dtype=float)
    spread = np.nanmax(P, axis=1) - np.nanmin(P, axis=1)
    distinctness = float(np.nanmean(spread > distinct_thresh))

    ch = cl.get("ch_scores", {})
    ks = sorted(int(k) for k in ch)
    chosen = int(cl["k"])
    ch_boundary = bool(ks and chosen == max(ks) and ch[chosen] == max(ch.values()))

    artifact = bool(distinctness < distinct_floor or ch_boundary)
    res = {"clone_cna_distinctness": distinctness, "chosen_k": chosen,
           "ch_boundary_k": ch_boundary,
           "clones_likely_artifact": artifact,
           "verdict": ("clones lack distinct arm-level CNA / k pinned to boundary "
                       "-> likely coverage/smoothing artifact, NOT real subclones" if artifact
                       else "clones show distinct CNA and CH chose k on merit")}
    adata.uns["spatial_omics_clone_diagnostics"] = res
    return res


def detect_channel_stripes(adata, key, *, ratio=3.0):
    """Detect grid-aligned (full-row / full-column) banding in a per-spot map.

    DBiT coverage is dominated by per-channel efficiency: a set of low-efficiency X- or
    Y-channels shows up as full rows/columns of a per-spot map (e.g. a region z-score).
    This is a microfluidic artifact, not spatial biology. We compare the variance of
    per-row (and per-column) means against the residual after removing those band means;
    a large ratio means the map is band-dominated.

    ``key`` is a column in ``adata.obs`` (e.g. ``'z_chr8'`` or ``'total_frags'``).
    Returns a dict with row/column band-dominance ratios and a ``stripe_artifact`` flag.
    """
    obs = adata.obs
    if key not in obs:
        raise KeyError(f"{key!r} not in adata.obs")
    x = obs["x_id"].values.astype(int); y = obs["y_id"].values.astype(int)
    v = obs[key].values.astype(float)
    nx, ny = x.max(), y.max()
    grid = np.full((nx, ny), np.nan)
    grid[x - 1, y - 1] = v
    finite = np.isfinite(grid)
    total_var = np.nanvar(grid)
    if total_var == 0:
        return {"row_dominance": 0.0, "col_dominance": 0.0, "stripe_artifact": False}
    row_means = np.nanmean(grid, axis=1, keepdims=True)
    col_means = np.nanmean(grid, axis=0, keepdims=True)
    row_resid = np.nanvar(grid - row_means)
    col_resid = np.nanvar(grid - col_means)
    row_dom = float((total_var - row_resid) / max(total_var, 1e-12))
    col_dom = float((total_var - col_resid) / max(total_var, 1e-12))
    # 'dominance' is the fraction of variance explained by pure row/column band means
    stripe = bool(max(row_dom, col_dom) > 0.5 and
                  (row_dom / max(col_dom, 1e-6) > ratio or col_dom / max(row_dom, 1e-6) > ratio))
    res = {"row_dominance": row_dom, "col_dominance": col_dom,
           "stripe_artifact": stripe,
           "axis": ("X-channels (rows)" if row_dom > col_dom else "Y-channels (columns)"),
           "verdict": ("grid-aligned banding -> likely microfluidic channel artifact"
                       if stripe else "no dominant channel banding")}
    adata.uns[f"spatial_omics_stripes_{key}"] = res
    return res
