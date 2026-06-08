"""Tools: the analysis core for sparse spatial DNA-seq.

    dual_smooth -> call_clones -> copy_number
                \-> permutation_significance

Three methods carried over from slide-DNA-seq (Zhao et al. Nature 2022) and the
stage3 pipeline, the parts that make sparse spatial genomics work:

* **dual_smooth** — PC-space + xy-space double smoothing of PC scores; denoises
  without blurring real clone boundaries (paper ``pc_scores_smo_both``).
* **permutation_significance** — empirical-null z-scores by reassigning fragments
  to spots, giving a spatially-resolved significance map (paper Fig 2c).
* **call_clones** — de-novo PCA + KMeans with Calinski–Harabasz model selection.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score
from sklearn.neighbors import NearestNeighbors

from .._constants import autosomes, chrom_order

__all__ = ["dual_smooth", "call_clones", "copy_number", "permutation_significance",
           "variance_decomposition"]


# ----------------------------------------------------------------------------- helpers
def _present(adata):
    return adata.obs["total_frags"].values > 0


def _qc_idx(adata):
    if "pass_qc" in adata.var:
        return np.where(adata.var["pass_qc"].values)[0]
    return np.arange(adata.n_vars)


def _genome_smooth(L, chrom_of_bin, win):
    """Median-filter along the genome within each chromosome (kept-bin order)."""
    out = L.copy()
    for c in pd.unique(chrom_of_bin):
        idx = np.where(chrom_of_bin == c)[0]
        if len(idx) >= win:
            out[:, idx] = median_filter(L[:, idx], size=(1, win), mode="nearest")
    return out


def _double_smooth_pcs(pcs, xy, k_pc=50, k_xy=50, k_join=10):
    """Smooth PC scores in PC-space, then average over xy neighbours (paper trick)."""
    n = len(pcs)
    kpc, kxy = min(k_pc, n), min(k_xy, n)
    _, idx_pc = NearestNeighbors(n_neighbors=kpc).fit(pcs).kneighbors(pcs)
    smo_pc = pcs[idx_pc].mean(axis=1)
    _, idx_xy = NearestNeighbors(n_neighbors=kxy).fit(xy).kneighbors(xy)
    kj = min(k_join, idx_xy.shape[1])
    return smo_pc[idx_xy[:, :kj]].mean(axis=1)


def _choose_k(X, kmin, kmax, seed):
    X = np.nan_to_num(np.asarray(X, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    kmax = min(kmax, X.shape[0] - 1)
    scores, best = {}, (kmin, -1.0, None)
    for k in range(kmin, kmax + 1):
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
        if len(np.unique(km.labels_)) < 2:
            continue  # degenerate solution — CH undefined
        s = calinski_harabasz_score(X, km.labels_)
        scores[k] = float(s)
        if s > best[1]:
            best = (k, s, km.labels_)
    if best[2] is None:  # fallback: force a 2-way split
        labels = KMeans(n_clusters=2, n_init=10, random_state=seed).fit_predict(X)
        return 2, labels, {2: 0.0}
    return best[0], best[2], scores


# ----------------------------------------------------------------------------- dual_smooth
def dual_smooth(adata, *, n_pcs=15, k_pc=50, k_xy=50, k_join=10,
                genome_win=3, regress_depth=True, seed=0):
    """log2-ratio -> depth-residualize -> genome-smooth -> PCA -> double-smooth PCs.

    Writes ``obsm['X_pca']``, ``obsm['X_pca_smoothed']`` and ``layers['log2_ratio']``
    (full-var, NaN outside QC bins). QC of the depth confound (corr of PC1 with
    log coverage) goes to ``uns['spatial_omics_dual_smooth']``.
    """
    pres = _present(adata)
    keep = _qc_idx(adata)
    rel = np.asarray(adata.layers["relative"])[:, keep]

    L = np.log2(np.clip(rel, 1e-3, None))
    L = np.nan_to_num(L, nan=0.0)
    L = L - np.median(L, axis=1, keepdims=True)  # remove per-spot ploidy/DC offset

    tot = adata.obs["total_frags"].values.astype(float)
    if regress_depth:
        cov = np.log10(tot + 1.0)
        covp = cov[pres] - cov[pres].mean()
        denom = float((covp ** 2).sum())
        if denom > 0:
            Lp = L[pres]
            beta = (Lp * covp[:, None]).sum(axis=0) / denom
            L[pres] = Lp - covp[:, None] * beta[None, :]

    chrom_of_keep = adata.var["chr"].values[keep]
    L = _genome_smooth(L, chrom_of_keep, genome_win)
    L = np.clip(L, -2.0, 2.0)

    # scatter back to full-var for downstream profile alignment
    full = np.full((adata.n_obs, adata.n_vars), np.nan, dtype=np.float32)
    full[:, keep] = L
    adata.layers["log2_ratio"] = full

    Lp = np.nan_to_num(L[pres], nan=0.0, posinf=0.0, neginf=0.0)
    npc = min(n_pcs, Lp.shape[1], max(Lp.shape[0] - 1, 1))
    pcs = PCA(n_components=npc, random_state=seed).fit_transform(Lp)
    pcs = np.nan_to_num(pcs, nan=0.0, posinf=0.0, neginf=0.0)
    xy = adata.obsm["spatial"][pres]
    pcs_smo = _double_smooth_pcs(pcs, xy, k_pc, k_xy, k_join)

    Xp = np.full((adata.n_obs, npc), np.nan, dtype=np.float32)
    Xs = np.full((adata.n_obs, npc), np.nan, dtype=np.float32)
    Xp[pres], Xs[pres] = pcs, pcs_smo
    adata.obsm["X_pca"] = Xp
    adata.obsm["X_pca_smoothed"] = Xs

    r_depth = float(np.corrcoef(pcs[:, 0], np.log10(tot[pres] + 1))[0, 1])
    adata.uns["spatial_omics_dual_smooth"] = {"n_pcs": int(npc), "depth_confound_r_pc1": r_depth}
    return adata


# ----------------------------------------------------------------------------- call_clones
def call_clones(adata, *, kmin=2, kmax=8, seed=0):
    """De-novo clones via KMeans on smoothed PCs (Calinski–Harabasz selects k).

    The flattest cluster (lowest median |log2-ratio|) is labelled ``Normal``;
    the rest are ``Clone1, Clone2, ...`` ordered by mean aberration. Writes
    ``obs['cluster']`` / ``obs['clone']`` and per-clone CNA profiles to
    ``uns['spatial_omics_clones']``.
    """
    if "X_pca_smoothed" not in adata.obsm:
        raise ValueError("run spatial_omics.tl.dual_smooth(adata) first")
    pres = _present(adata)
    X = adata.obsm["X_pca_smoothed"][pres]
    k, lab_p, scores = _choose_k(X, kmin, kmax, seed)

    labels = np.full(adata.n_obs, -1)
    labels[pres] = lab_p

    L = adata.layers["log2_ratio"]
    Lp = np.nan_to_num(L[pres], nan=0.0)
    flat = {c: np.median(np.abs(Lp[lab_p == c])) for c in range(k)}
    normal_c = min(flat, key=flat.get)

    clone_mean = np.vstack([np.nanmean(L[labels == c], axis=0) for c in range(k)])
    cna = clone_mean - clone_mean[normal_c][None, :]
    cna = cna - np.nanmedian(cna, axis=1, keepdims=True)

    aber = {c: np.nanmean(np.abs(cna[c])) for c in range(k)}
    tumor = sorted([c for c in range(k) if c != normal_c], key=lambda c: -aber[c])
    name = {normal_c: "Normal"}
    for i, c in enumerate(tumor, 1):
        name[c] = f"Clone{i}"

    adata.obs["cluster"] = labels
    adata.obs["clone"] = pd.Categorical(
        [name.get(c, "NA") if c >= 0 else "empty" for c in labels])

    profiles = pd.DataFrame(cna.T, index=adata.var_names,
                            columns=[name[c] for c in range(k)])
    adata.uns["spatial_omics_clones"] = {
        "k": int(k), "normal_cluster": int(normal_c),
        "names": {int(c): name[c] for c in range(k)},
        "ch_scores": scores, "profiles": profiles,
        "sizes": {name[c]: int((labels == c).sum()) for c in range(k)},
    }
    return adata


# ----------------------------------------------------------------------------- copy_number
def copy_number(adata, *, k_spatial=49, genome_win=3, cna_thresh=0.5, normal_mask=None):
    """Per-spot integer-scaled copy number via spatial-neighbourhood densification.

    Sums corrected counts over each spot's ``k_spatial`` nearest neighbours, library-
    normalizes, divides by a per-bin pseudo-normal reference, anchors the autosomal
    mode to 2, genome-smooths and clips to [0, 6]. Writes ``layers['copy_number']``
    (full-var) and ``obs['cnv_burden']``. Mirrors stage3 per-spot CNV.
    """
    pres = _present(adata)
    keep = _qc_idx(adata)
    layer = "corrected" if "corrected" in adata.layers else "counts"
    C = np.asarray(adata.layers[layer], dtype=float)[:, keep]

    xy = adata.obsm["spatial"]
    kk = min(k_spatial, int(pres.sum()))
    nn = NearestNeighbors(n_neighbors=kk).fit(xy[pres]).kneighbors(xy[pres], return_distance=False)
    Cp = C[pres]
    Cs_p = Cp[nn].sum(axis=1)  # neighbourhood SUM (densify)
    Cs = np.zeros_like(C)
    Cs[pres] = Cs_p

    t2 = Cs.sum(1)
    med = np.median(t2[pres]) if pres.any() else 1.0
    libn = np.zeros_like(Cs)
    libn[pres] = Cs[pres] / t2[pres, None] * med

    ref_spots = pres
    if normal_mask is not None:
        nm = pres & np.asarray(normal_mask, bool)
        if nm.sum() >= 20:
            ref_spots = nm
    ref = np.median(libn[ref_spots], axis=0)
    ref[ref <= 0] = np.nan
    rel = np.full_like(libn, np.nan)
    rel[pres] = libn[pres] / ref

    chrom_keep = adata.var["chr"].values[keep]
    auto = set(autosomes())
    is_auto = np.array([c in auto for c in chrom_keep])
    base = np.nanmedian(rel[:, is_auto], axis=1, keepdims=True)
    base[base <= 0] = np.nan
    CN = rel / base * 2.0
    CN = _genome_smooth(np.nan_to_num(CN, nan=2.0), chrom_keep, genome_win)
    CN = np.clip(CN, 0, 6)

    full = np.full((adata.n_obs, adata.n_vars), np.nan, dtype=np.float32)
    full[:, keep] = CN
    adata.layers["copy_number"] = full
    burden = np.full(adata.n_obs, np.nan)
    burden[pres] = (np.abs(CN[pres] - 2) > cna_thresh).mean(1)
    adata.obs["cnv_burden"] = burden
    return adata


# ----------------------------------------------------------------------------- significance
def permutation_significance(adata, *, regions=("chr8", "chr17", "chr18"),
                             n_perm=100, k_xy=10, seed=0, normal_label="Normal"):
    """Spatially-resolved significance per genomic region (paper Fig 2c).

    For each region, reassign all its fragments to spots by a multinomial null
    (p ~ each spot's fragment share), preserving the region total, over ``n_perm``
    draws; xy-smooth (sum over ``k_xy`` neighbours); z=(obs-mu)/sd; subtract the
    Normal-cluster mean. Writes ``obs['z_<region>']``.
    """
    pres = _present(adata)
    idx = np.where(pres)[0]
    counts = np.asarray(adata.layers["counts"])
    tot = adata.obs["total_frags"].values.astype(float)
    xy = adata.obsm["spatial"][pres]

    if "clone" in adata.obs:
        normal_mask = (adata.obs["clone"].values == normal_label)
    else:
        normal_mask = np.zeros(adata.n_obs, bool)

    kk = min(k_xy, len(idx))
    nn = NearestNeighbors(n_neighbors=kk).fit(xy).kneighbors(xy, return_distance=False)
    rng = np.random.default_rng(seed)
    chrom = adata.var["chr"].values

    for region in regions:
        sel = (chrom == region) if region.startswith("chr") else _parse_region(adata, region)
        if not np.any(sel):
            continue
        rc = counts[:, sel].sum(axis=1)
        obs = rc[idx].astype(float)
        tt = tot[idx]
        M = int(round(obs.sum()))
        p = tt / tt.sum() if tt.sum() > 0 else np.ones(len(tt)) / len(tt)
        perm = rng.multinomial(M, p, size=n_perm).astype(float)
        obs_s = obs[nn].sum(axis=1)
        perm_s = perm[:, nn].sum(axis=2)
        mu, sd = perm_s.mean(0), perm_s.std(0) + 1e-9
        z = (obs_s - mu) / sd
        nm = normal_mask[idx]
        if nm.any():
            z = z - np.nanmean(z[nm])
        col = np.full(adata.n_obs, np.nan)
        col[idx] = z
        adata.obs[f"z_{region}"] = col
    return adata


def _parse_region(adata, region):
    """Accept 'chr8:120-145' (Mb) region strings."""
    chrom, rng = region.split(":")
    lo, hi = (float(v) * 1e6 for v in rng.split("-"))
    v = adata.var
    return (v["chr"].values == chrom) & (v["start"].values >= lo) & (v["start"].values < hi)


# ----------------------------------------------------------------------------- variance decomposition
def _align_spots(dna, rna):
    """Match spots between DNA and RNA AnnData by (x_id, y_id). Returns index arrays."""
    dkey = {(int(x), int(y)): i for i, (x, y) in enumerate(zip(dna.obs.x_id, dna.obs.y_id))}
    di, ri = [], []
    for j, (x, y) in enumerate(zip(rna.obs.x_id, rna.obs.y_id)):
        i = dkey.get((int(x), int(y)))
        if i is not None:
            di.append(i)
            ri.append(j)
    return np.array(di, dtype=int), np.array(ri, dtype=int)


def _rss(X, Y):
    """Residual sum of squares per column of Y for OLS fit Y ~ X (vectorized)."""
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    return (resid ** 2).sum(axis=0)


def variance_decomposition(dna, rna, *, clone_key="clone", celltype_key="cell_type",
                           tumor_types=("Tumor",), immune_types=("Immune",),
                           k_density=50, layer=None, normal_label="Normal"):
    """Partition each gene's spatial expression variance into genetic vs microenvironment.

    The multi-modal core (paper Fig 4). For spots shared by the DNA and RNA assays:

    * **genetic (intrinsic)** = variance explained by DNA **subclone identity**
      (``dna.obs[clone_key]``);
    * **microenvironment (extrinsic)** = variance explained by local **tumour density**
      and **immune density** (fraction of each spot's ``k_density`` RNA neighbours that
      are tumour / immune cells, from ``rna.obs[celltype_key]``).

    For every gene an OLS model ``expr ~ clone + tumor_density + immune_density`` is fit;
    each factor's contribution is its extra sum of squares (drop in RSS vs the model
    without it), normalized by total SS. Writes per-gene fractions to ``rna.var`` and a
    genome-wide summary to ``rna.uns['spatial_omics_vardecomp']``.

    Returns
    -------
    The ``rna`` AnnData, annotated in place.
    """
    di, ri = _align_spots(dna, rna)
    if len(di) < 20:
        raise ValueError(f"only {len(di)} spots shared between DNA and RNA assays")

    # expression matrix (shared spots x genes)
    Y = rna.layers[layer] if layer else rna.X
    Y = np.asarray(Y[ri], dtype=float)
    Y = Y - Y.mean(axis=0, keepdims=True)
    tss = (Y ** 2).sum(axis=0)
    tss[tss <= 0] = np.nan

    n = len(ri)
    intercept = np.ones((n, 1))

    # genetic factor: clone one-hot (drop-first), from DNA on shared spots
    clones = pd.Series(dna.obs[clone_key].astype(str).values[di])
    clone_dummies = pd.get_dummies(clones, drop_first=True).values.astype(float)
    if clone_dummies.shape[1] == 0:
        clone_dummies = np.zeros((n, 0))

    # microenvironment factors: tumour / immune density from RNA cell types
    ct = rna.obs[celltype_key].astype(str).values[ri]
    xy = rna.obsm["spatial"][ri]
    kk = min(k_density, n)
    nn = NearestNeighbors(n_neighbors=kk).fit(xy).kneighbors(xy, return_distance=False)
    is_tumor = np.isin(ct, tumor_types).astype(float)
    is_immune = np.isin(ct, immune_types).astype(float)
    tumor_density = is_tumor[nn].mean(axis=1)
    immune_density = is_immune[nn].mean(axis=1)
    dens = np.column_stack([tumor_density - tumor_density.mean(),
                            immune_density - immune_density.mean()])

    X_full = np.hstack([intercept, clone_dummies, dens])
    X_no_clone = np.hstack([intercept, dens])
    X_no_dens = np.hstack([intercept, clone_dummies])

    rss_full = _rss(X_full, Y)
    ss_genetic = np.clip(_rss(X_no_clone, Y) - rss_full, 0, None)
    ss_micro = np.clip(_rss(X_no_dens, Y) - rss_full, 0, None)

    frac_gen = ss_genetic / tss
    frac_mic = ss_micro / tss
    frac_res = rss_full / tss

    rna.var["vd_genetic"] = frac_gen
    rna.var["vd_microenv"] = frac_mic
    rna.var["vd_residual"] = frac_res
    rna.var["vd_dominant"] = np.where(frac_gen >= frac_mic, "genetic", "microenv")

    rna.uns["spatial_omics_vardecomp"] = {
        "n_spots_shared": int(n),
        "mean_frac_genetic": float(np.nanmean(frac_gen)),
        "mean_frac_microenv": float(np.nanmean(frac_mic)),
        "mean_frac_residual": float(np.nanmean(frac_res)),
        "n_genes_genetic_dominant": int(np.nansum(frac_gen > frac_mic)),
        "n_genes_microenv_dominant": int(np.nansum(frac_mic > frac_gen)),
    }
    return rna

# rigor layer — confound controls / artifact guards (the honest core)
from .rigor import (  # noqa: E402
    morans_i,
    spatial_heterogeneity,
    clone_diagnostics,
    detect_channel_stripes,
)

# H&E morphology — independent purity/cellularity (optional [he] extra)
from .he import he_purity  # noqa: E402
