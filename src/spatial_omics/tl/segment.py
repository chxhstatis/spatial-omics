"""CNA segmentation — piecewise-constant copy-number segments from bin-level values.

Bin-level copy number is noisy. Segmentation partitions each chromosome into runs of
constant copy number (the basis for integer-CN calling, breakpoints, and focal-event
detection). We use recursive **binary segmentation** with a BIC-style penalty (the core
of CBS, dependency-free: numpy only). Segments never cross chromosome boundaries.

Guard: reports variance-explained and a ``flat_genome`` flag — on a low-purity / event-free
profile, segmentation honestly finds ~one segment per chromosome (no CNA), not spurious ones.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .._constants import chrom_order


def _best_split(y, min_size):
    """index k (and SSE reduction) of the best single split of 1D y, or None."""
    n = len(y)
    if n < 2 * min_size:
        return None
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y * y)])
    ks = np.arange(min_size, n - min_size + 1)
    sumL = cs[ks]; sumR = cs[n] - sumL
    sseL = cs2[ks] - sumL ** 2 / ks
    sseR = (cs2[n] - cs2[ks]) - sumR ** 2 / (n - ks)
    total = cs2[n] - cs[n] ** 2 / n
    red = total - sseL - sseR
    i = int(np.argmax(red))
    return int(ks[i]), float(red[i])


def _binseg(y, penalty, min_size):
    """recursive binary segmentation -> sorted breakpoint indices within y."""
    cps = []
    stack = [(0, len(y))]
    while stack:
        a, b = stack.pop()
        r = _best_split(y[a:b], min_size)
        if r is None:
            continue
        k, red = r
        if red > penalty:
            cps.append(a + k)
            stack.append((a, a + k))
            stack.append((a + k, b))
    return sorted(cps)


def segment(adata, *, layer="copy_number", profile=None, penalty_scale=2.0, min_size=3):
    """Segment the (pseudobulk) copy-number profile into piecewise-constant runs.

    Parameters
    ----------
    layer : str
        Per-spot layer to pseudobulk (median over present spots) when ``profile`` is None.
    profile : array (n_vars,), optional
        Segment this 1-D profile directly (NaN allowed) instead of the pseudobulk — e.g. a
        single clone's profile. Order must match ``adata.var``.
    penalty_scale : float
        Multiplies the BIC-style penalty ``sigma^2 * log(N)`` for accepting a split. Higher
        = fewer, larger segments. Noise ``sigma`` is estimated robustly from bin-to-bin diffs.
    min_size : int
        Minimum bins per segment.

    Writes ``adata.var['segment']`` (genome-wide segment id; -1 where unsegmented),
    ``adata.var['seg_mean']`` (segmented value broadcast to bins), and
    ``adata.uns['spatial_omics_segments']`` = {table, n_segments, var_explained, flat_genome}.
    Returns the segments DataFrame (chr, start, end, n_bins, mean).
    """
    var = adata.var
    if profile is None:
        M = np.asarray(adata.layers[layer], dtype=float)
        tot = np.asarray(adata.layers["counts"]).sum(1)
        prof = np.nanmedian(M[tot > 0], axis=0)
    else:
        prof = np.asarray(profile, dtype=float)
    if prof.shape[0] != adata.n_vars:
        raise ValueError("profile length must equal adata.n_vars")

    chrom = var["chr"].astype(str).values
    start = var["start"].values
    end = var["end"].values if "end" in var else start

    seg_id = np.full(adata.n_vars, -1, dtype=int)
    seg_mean = np.full(adata.n_vars, np.nan)
    rows = []
    sid = 0
    for c in chrom_order():
        idx = np.where((chrom == c) & np.isfinite(prof))[0]
        if len(idx) < min_size:
            continue
        y = prof[idx]
        d = np.diff(y)
        sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) if len(d) else 0.0
        penalty = penalty_scale * max(sigma, 1e-6) ** 2 * np.log(max(len(y), 2))
        cps = _binseg(y, penalty, min_size)
        bounds = [0] + cps + [len(y)]
        for a, b in zip(bounds[:-1], bounds[1:]):
            members = idx[a:b]
            m = float(np.mean(y[a:b]))
            seg_id[members] = sid
            seg_mean[members] = m
            rows.append({"segment": sid, "chr": c, "start": int(start[members[0]]),
                         "end": int(end[members[-1]]), "n_bins": int(b - a), "mean": round(m, 4)})
            sid += 1

    fin = np.isfinite(prof) & np.isfinite(seg_mean)
    total_var = float(np.var(prof[fin])) if fin.any() else 0.0
    resid_var = float(np.var(prof[fin] - seg_mean[fin])) if fin.any() else 0.0
    var_explained = float(1 - resid_var / total_var) if total_var > 0 else 0.0
    n_chrom = len(set(chrom[np.isfinite(prof)]))
    flat = bool(sid <= n_chrom)

    adata.var["segment"] = seg_id
    adata.var["seg_mean"] = seg_mean.astype(np.float32)
    table = pd.DataFrame(rows)
    adata.uns["spatial_omics_segments"] = {
        "table": table, "n_segments": int(sid),
        "var_explained": round(var_explained, 3), "flat_genome": flat,
        "verdict": ("near-flat profile: ~one segment per chromosome (no CNA segments) — "
                    "consistent with low purity / no events" if flat else
                    "%d copy-number segments called" % sid)}
    return table
