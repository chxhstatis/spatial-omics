"""Plotting: the figures that sell the method.

Every function returns a matplotlib ``Figure`` so it composes in notebooks and the
docs site. Spatial maps assume a DBiT grid via ``obsm['spatial']``.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg") if matplotlib.get_backend().lower() == "agg" else None
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from .._constants import chrom_order

__all__ = ["spatial_clones", "spatial_copy_number", "significance", "clone_profiles",
           "variance_decomposition"]


def _grid(adata):
    return int(adata.uns.get("spatial_omics", {}).get("grid", 50))


def _to_grid(adata, values):
    g = _grid(adata)
    arr = np.full((g, g), np.nan)
    xs = adata.obs["x_id"].values.astype(int) - 1
    ys = adata.obs["y_id"].values.astype(int) - 1
    arr[xs, ys] = values
    return arr


def spatial_clones(adata, *, key="clone", ax=None, title=None):
    """Paint the de-novo clone label on the grid."""
    cats = list(adata.obs[key].cat.categories) if hasattr(adata.obs[key], "cat") \
        else sorted(set(adata.obs[key]))
    code = {c: i for i, c in enumerate(cats)}
    vals = np.array([code[c] for c in adata.obs[key]], dtype=float)
    grid = _to_grid(adata, vals)
    fig = ax.figure if ax else plt.figure(figsize=(6, 5))
    ax = ax or plt.gca()
    cmap = plt.get_cmap("tab10", max(len(cats), 2))
    im = ax.imshow(grid, origin="lower", cmap=cmap, vmin=-0.5, vmax=len(cats) - 0.5,
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, ticks=range(len(cats)))
    cbar.ax.set_yticklabels(cats)
    ax.set_title(title or "De-novo clones"); ax.set_xlabel("Y"); ax.set_ylabel("X")
    return fig


def spatial_copy_number(adata, *, chrom="chr8", ax=None, title=None):
    """Mean copy number of a chromosome painted on the grid."""
    keep = adata.var["chr"].values == chrom
    cn = np.nanmean(np.asarray(adata.layers["copy_number"])[:, keep], axis=1)
    grid = _to_grid(adata, cn)
    fig = ax.figure if ax else plt.figure(figsize=(6, 5))
    ax = ax or plt.gca()
    im = ax.imshow(grid, origin="lower", cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=1, vcenter=2, vmax=3), interpolation="nearest")
    fig.colorbar(im, ax=ax, label="copy number")
    ax.set_title(title or f"Copy number — {chrom}"); ax.set_xlabel("Y"); ax.set_ylabel("X")
    return fig


def significance(adata, *, region="chr8", ax=None, title=None):
    """Signed permutation z-score map for a region (red=gain, blue=loss)."""
    col = f"z_{region}"
    if col not in adata.obs:
        raise ValueError(f"run tl.permutation_significance(adata, regions=['{region}']) first")
    grid = _to_grid(adata, adata.obs[col].values)
    vmax = np.nanmax(np.abs(grid)) or 1.0
    fig = ax.figure if ax else plt.figure(figsize=(6, 5))
    ax = ax or plt.gca()
    im = ax.imshow(grid, origin="lower", cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), interpolation="nearest")
    fig.colorbar(im, ax=ax, label="z (vs permutation null)")
    ax.set_title(title or f"Spatial significance — {region}"); ax.set_xlabel("Y"); ax.set_ylabel("X")
    return fig


def clone_profiles(adata, *, figsize=(13, None)):
    """Genome-wide CNA profile per clone (from uns['spatial_omics_clones'])."""
    res = adata.uns["spatial_omics_clones"]
    prof = res["profiles"]
    chrom = adata.var["chr"].values
    bnds, centers, order, pos = [], [], [], 0
    for c in chrom_order():
        n = int((chrom == c).sum())
        if n == 0:
            continue
        bnds.append(pos); centers.append(pos + n / 2); order.append(c.replace("chr", "")); pos += n
    cols = [c for c in prof.columns if c != "Normal"] or list(prof.columns)
    h = figsize[1] or 1.6 * len(cols)
    fig, axes = plt.subplots(len(cols), 1, figsize=(figsize[0], h), squeeze=False)
    for ax, name in zip(axes[:, 0], cols):
        ax.plot(prof[name].values, lw=0.5, color="#b2182b")
        ax.axhline(0, color="grey", lw=0.5)
        for b in bnds:
            ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.4)
        ax.set_xticks(centers); ax.set_xticklabels(order, fontsize=7)
        ax.set_ylabel(f"{name}\nlog2 ratio", fontsize=8)
    axes[0, 0].set_title("Clone copy-number profiles (vs Normal)")
    fig.tight_layout()
    return fig


def variance_decomposition(rna, *, n_top=40, ax=None):
    """Stacked bars of genetic vs microenvironment vs residual variance per gene.

    Shows the ``n_top`` genes with the most explained (genetic+microenv) variance,
    the figure that makes the intrinsic/extrinsic split legible (paper Fig 4).
    """
    v = rna.var
    expl = (v["vd_genetic"].fillna(0) + v["vd_microenv"].fillna(0)).values
    order = np.argsort(-expl)[:n_top]
    g = v["vd_genetic"].values[order]
    m = v["vd_microenv"].values[order]
    res = v["vd_residual"].values[order]
    names = v.index.values[order]
    fig = ax.figure if ax else plt.figure(figsize=(max(8, n_top * 0.22), 4))
    ax = ax or plt.gca()
    x = np.arange(len(order))
    ax.bar(x, g, label="genetic (clone)", color="#b2182b")
    ax.bar(x, m, bottom=g, label="microenvironment", color="#2166ac")
    ax.bar(x, res, bottom=g + m, label="residual", color="#cccccc")
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel("fraction of variance"); ax.set_ylim(0, 1)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Expression variance: intrinsic (genetic) vs extrinsic (microenvironment)")
    fig.tight_layout()
    return fig
