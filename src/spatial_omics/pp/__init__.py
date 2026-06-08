"""Preprocessing: bias correction & normalization for sparse spatial DNA-seq.

Coverage is NOT copy number — between them sit GC, mappability and depth biases.
This module reproduces the stage3 correction order:

    add_reference -> bin_qc -> correct_bias (GC + mappability) -> normalize

after which ``tl`` can smooth, call clones, and derive copy number.
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from .._bins import bin_offsets
from .._constants import autosomes

__all__ = ["add_reference", "bin_qc", "correct_bias", "normalize",
           "load_reference_tracks", "default_reference_dir",
           "normal_anchor", "pick_normal_spots"]


def default_reference_dir() -> str:
    """Path to the public bias tracks bundled inside the package (hg38, 1 Mb)."""
    import os
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "_ref")


def load_reference_tracks(adata, ref_dir: str = None, *, bin_size=1_000_000, genome="hg38"):
    """Attach bias tracks (GC / mappability / blacklist) to ``adata.var``.

    ``ref_dir`` defaults to the public tracks bundled with the package
    (``default_reference_dir()``), so real-data runs work out of the box. Pass a
    custom dir to use your own ``bin_gc_1mb.tsv`` / ``k100.umap.bed.gz`` /
    ``hg38-blacklist.v2.bed.gz``.
    """
    import os

    if ref_dir is None:
        ref_dir = default_reference_dir()
    off, n_bins = bin_offsets(bin_size, genome)
    chroms = set(adata.var["chr"])

    gcf = pd.read_csv(os.path.join(ref_dir, "bin_gc_1mb.tsv"), sep="\t")
    gmap = dict(zip(gcf["chr"] + ":" + (gcf["start"] // bin_size).astype(str), gcf["gc"]))
    key = adata.var["chr"] + ":" + (adata.var["start"] // bin_size).astype(str)
    gc = key.map(gmap).fillna(0.0).values

    def _bed_frac(path, skip_track):
        cov = np.zeros(n_bins)
        op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        with op as f:
            for line in f:
                if skip_track and line.startswith("track"):
                    continue
                p = line.split("\t")
                if len(p) < 3 or p[0] not in chroms:
                    continue
                c, s, e = p[0], int(p[1]), int(p[2])
                for b in range(s // bin_size, (e - 1) // bin_size + 1):
                    cov[off[c] + b] += min(e, (b + 1) * bin_size) - max(s, b * bin_size)
        return cov / bin_size

    mapp = _bed_frac(os.path.join(ref_dir, "k100.umap.bed.gz"), True)
    blk = _bed_frac(os.path.join(ref_dir, "hg38-blacklist.v2.bed.gz"), False)
    return add_reference(adata, gc=gc, mappability=mapp, blacklist=blk)


def add_reference(adata, *, gc=None, mappability=None, blacklist=None):
    """Attach per-bin bias covariates to ``adata.var`` (gc / mappability / blacklist)."""
    if gc is not None:
        adata.var["gc"] = np.asarray(gc, dtype=float)
    if mappability is not None:
        adata.var["mappability"] = np.asarray(mappability, dtype=float)
    if blacklist is not None:
        adata.var["blacklist"] = np.asarray(blacklist, dtype=float)
    return adata


def bin_qc(adata, *, gc_min=0.35, map_min=0.5, blacklist_max=0.3, min_bin_frac=0.2):
    """Flag bins that pass QC (-> ``adata.var['pass_qc']``).

    Mirrors stage3: drop low-GC, low-mappability, blacklisted, and under-covered bins.
    """
    counts = adata.layers["counts"]
    pres = counts.sum(1) > 0
    bin_mean = counts[pres].mean(0)
    med = np.median(bin_mean[bin_mean > 0]) if (bin_mean > 0).any() else 0.0
    keep = bin_mean > (min_bin_frac * med)
    if "gc" in adata.var:
        keep &= adata.var["gc"].values > gc_min
    if "mappability" in adata.var:
        keep &= adata.var["mappability"].values > map_min
    if "blacklist" in adata.var:
        keep &= adata.var["blacklist"].values < blacklist_max
    adata.var["pass_qc"] = keep
    return adata


def _qcorrect_rows(C, covar, n_q=20):
    """Per-spot flatten of counts against a bin-level covariate (GC or mappability).

    Bins are split into ``n_q`` equal-count quantile groups by ``covar``; each spot's
    counts are divided by the spot's mean count in that group, then rescaled to the
    spot's overall mean. Ported from stage3 ``qcorrect_rows``.
    """
    covar = np.asarray(covar, dtype=float)
    ranks = np.argsort(np.argsort(covar))
    grp = np.minimum(ranks * n_q // len(covar), n_q - 1)
    oh = np.zeros((len(covar), n_q))
    oh[np.arange(len(covar)), grp] = 1.0
    gcount = oh.sum(0)
    gsum = C @ oh
    expected = np.divide(gsum, gcount[None, :], out=np.zeros_like(gsum), where=gcount[None, :] > 0)
    exp_bin = expected[:, grp]
    rowmean = C.mean(1, keepdims=True)
    return np.divide(C, exp_bin, out=np.zeros_like(C), where=exp_bin > 0) * rowmean


def correct_bias(adata, *, n_quantiles=20):
    """GC + mappability quantile correction -> ``adata.layers['corrected']``."""
    C = np.asarray(adata.layers["counts"], dtype=float)
    if "gc" in adata.var:
        C = _qcorrect_rows(C, adata.var["gc"].values, n_quantiles)
    if "mappability" in adata.var:
        C = _qcorrect_rows(C, adata.var["mappability"].values, n_quantiles)
    adata.layers["corrected"] = C.astype(np.float32)
    return adata


def _spatial_densify(C, coords, present, k):
    """Neighbourhood-MEAN smoothing on the spatial grid (densifies sparse counts).

    Sparse spatial DNA-seq has ~1 fragment / bin / spot, so the cross-spot median is
    0 for almost every bin. Averaging each spot over its k nearest spatial neighbours
    densifies the matrix before the reference is built — the stage3 ``spatial_smooth``
    step. Operates on present spots only; absent spots stay 0.
    """
    out = np.zeros_like(C)
    idx = np.where(present)[0]
    if len(idx) == 0:
        return out
    kk = min(k, len(idx))
    nn = NearestNeighbors(n_neighbors=kk).fit(coords[idx]).kneighbors(coords[idx], return_distance=False)
    out[idx] = C[idx][nn].mean(axis=1)
    return out


def normalize(adata, *, normal_mask=None, spatial_k=49):
    """Spatially densify, library-size normalize, then divide by a pseudo-normal ref.

    Each spot is first averaged over its ``spatial_k`` nearest spatial neighbours
    (without this, the cross-spot median is 0 at typical sparse coverage). Reference =
    median over reference spots (``normal_mask`` if given, else all spots) on QC bins.
    Result -> ``adata.layers['relative']`` (NaN where a spot has no fragments).
    Ported from stage3 ``process``. Set ``spatial_k=1`` to disable densification.
    """
    layer = "corrected" if "corrected" in adata.layers else "counts"
    C = np.asarray(adata.layers[layer], dtype=float)
    tot0 = np.asarray(adata.layers["counts"]).sum(1)
    pres = tot0 > 0
    if spatial_k and spatial_k > 1:
        C = _spatial_densify(C, adata.obsm["spatial"], pres, spatial_k)
    tot = C.sum(1)
    med = np.median(tot[pres]) if pres.any() else 1.0
    libn = np.zeros_like(C)
    libn[pres] = C[pres] / tot[pres, None] * med

    pass_qc = adata.var["pass_qc"].values if "pass_qc" in adata.var else np.ones(adata.n_vars, bool)
    ref_spots = pres
    if normal_mask is not None:
        nm = pres & np.asarray(normal_mask, bool)
        if nm.sum() >= 20:
            ref_spots = nm
    ref = np.full(adata.n_vars, np.nan)
    ref[pass_qc] = np.median(libn[np.ix_(ref_spots, pass_qc)], axis=0)
    ref[ref <= 0] = np.nan

    rel = np.full_like(libn, np.nan)
    rel[pres] = libn[pres] / ref
    adata.layers["relative"] = rel.astype(np.float32)
    return adata


# ---------------------------------------------------------------- normal anchoring + guard
def _pseudobulk_log2_sd(adata):
    """SD across QC bins of the per-bin pseudobulk log2(relative) — the bulk CNA amplitude."""
    rel = np.asarray(adata.layers["relative"], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2 = np.log2(np.clip(rel, 1e-3, None))
    pb = np.nanmedian(log2, axis=0)
    pass_qc = adata.var["pass_qc"].values if "pass_qc" in adata.var else np.ones(adata.n_vars, bool)
    return float(np.nanstd(pb[pass_qc]))


def normal_anchor(adata, normal_mask, *, collapse_frac=0.4, spatial_k=49):
    """Re-normalize against a NORMAL/reference spot set, with a signal-collapse guard.

    Divides the per-spot profiles by a pseudo-normal built from ``normal_mask`` spots only
    (the paper's ÷bulk-WGS analogue). Crucially, it first measures whether anchoring to these
    spots *flattens* the copy-number signal: if the bulk CNA amplitude (pseudobulk log2 SD)
    drops below ``collapse_frac`` of its un-anchored value, the reference shares the tumour's
    own coverage pattern — there is NO genuine internal normal, and anchoring would erase real
    CNA. In that case use an EXTERNAL reference (control / Panel-of-Normals / RNA-defined stroma).

    The anchored ``layers['relative']`` is written regardless, but the diagnostic
    (``adata.uns['spatial_omics_normal_anchor']`` with ``signal_collapsed``) tells you whether
    to trust it. Returns ``adata``.
    """
    normalize(adata, normal_mask=None, spatial_k=spatial_k)
    sd_before = _pseudobulk_log2_sd(adata)
    normalize(adata, normal_mask=normal_mask, spatial_k=spatial_k)
    sd_after = _pseudobulk_log2_sd(adata)
    n_ref = int((np.asarray(normal_mask, bool) & (np.asarray(adata.layers["counts"]).sum(1) > 0)).sum())
    collapsed = bool(sd_before > 0 and sd_after < collapse_frac * sd_before)
    res = {"n_reference_spots": n_ref, "bulk_sd_before": sd_before, "bulk_sd_after": sd_after,
           "signal_collapsed": collapsed,
           "verdict": ("SIGNAL COLLAPSE: reference shares the tumour pattern -> no genuine internal "
                       "normal; use an EXTERNAL reference (control / RNA stroma). Anchored result not "
                       "trustworthy." if collapsed else
                       "normal-anchored OK (real CNA retained)")}
    adata.uns["spatial_omics_normal_anchor"] = res
    return adata


def pick_normal_spots(adata, *, quantile=0.30, min_frags=200, burden_key="cnv_burden"):
    """Data-driven candidate normal spots = the FLATTEST (lowest CNV-burden) spots.

    Requires ``tl.copy_number`` first (for ``cnv_burden``). Returns a boolean mask over spots.

    .. warning::
       "Flattest" is relative to the (possibly tumour-contaminated) cross-spot reference. In a
       uniformly low-purity section there is no true internal normal, and this set will share
       the tumour's pattern — anchoring to it erases the signal. Always feed the result through
       :func:`normal_anchor`, whose ``signal_collapsed`` guard detects exactly that failure.
       The reliable normal reference is EXTERNAL (control / RNA stroma).
    """
    if burden_key not in adata.obs:
        raise KeyError(f"{burden_key!r} not in adata.obs — run tl.copy_number first")
    burden = adata.obs[burden_key].values.astype(float)
    frags = adata.obs["total_frags"].values.astype(float)
    elig = (frags >= min_frags) & np.isfinite(burden)
    if elig.sum() < 20:
        elig = np.isfinite(burden) & (frags >= np.percentile(frags, 80))
    thr = np.quantile(burden[elig], quantile)
    return elig & (burden <= thr)
