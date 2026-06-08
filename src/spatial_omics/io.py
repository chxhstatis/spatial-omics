"""I/O: turn spatial_omics_pipeline outputs into a standard AnnData object.

The spatial_omics data object is a plain :class:`anndata.AnnData`:

    .X                     spot x bin counts (raw fragments)         [layer 'counts']
    .obs   (n_spots)       x_id, y_id, total_frags  (+ clone, burden after analysis)
    .var   (n_bins)        chr, start, end (+ gc, mappability, blacklist after pp)
    .obsm['spatial']       (n_spots, 2) physical xy coordinates on the grid
    .layers                'counts' / 'corrected' / 'copy_number' (filled by pp/tl)
    .obsp                  'spatial_knn' / 'pca_knn' (filled by tl.dual_smooth)
    .uns['spatial_omics']       sample id, grid, bin_size, genome, param/version provenance

Using AnnData means the whole scanpy / squidpy ecosystem works out of the box.
"""
from __future__ import annotations

import subprocess

import anndata as ad
import numpy as np
import pandas as pd

from ._bins import bin_offsets, build_bins
from ._constants import GENOMES, chrom_order

GRID = 50  # DBiT 50x50


def _read_fragments(frag_path: str) -> pd.DataFrame:
    cols = ["chr", "start", "end", "cb", "mapq"]
    if frag_path.endswith(".gz"):
        proc = subprocess.Popen(["gzip", "-dc", frag_path], stdout=subprocess.PIPE)
        fr = pd.read_csv(proc.stdout, sep="\t", header=None, names=cols, on_bad_lines="skip")
        proc.stdout.close()
        proc.wait()
    else:
        fr = pd.read_csv(frag_path, sep="\t", header=None, names=cols, on_bad_lines="skip")
    return fr


def from_pipeline(
    fragments: str,
    matched_barcodes: str,
    *,
    sample: str = "sample",
    bin_size: int = 1_000_000,
    genome: str = "hg38",
    grid: int = GRID,
) -> ad.AnnData:
    """Build an AnnData from Stage-2 outputs.

    Parameters
    ----------
    fragments
        ``top50x50_fragments_with_cb.tsv(.gz)`` — BED-like fragments with cell barcode.
    matched_barcodes
        ``matched_spot_barcodes.tsv`` with columns ``barcode_obs, x_id, y_id, full_status``.
    sample
        Sample name stored in ``.uns['spatial_omics']``.
    """
    sizes = GENOMES[genome]
    chroms = chrom_order(genome)
    off, n_bins = bin_offsets(bin_size, genome)
    bins = build_bins(bin_size, genome)

    m = pd.read_csv(matched_barcodes, sep="\t",
                    usecols=["barcode_obs", "x_id", "y_id", "full_status"])
    m = m[m["full_status"] == "mapped"].dropna(subset=["x_id", "y_id"])
    obs2x = dict(zip(m.barcode_obs, m.x_id.astype(int)))
    obs2y = dict(zip(m.barcode_obs, m.y_id.astype(int)))

    fr = _read_fragments(fragments)
    fr = fr[fr["chr"].isin(chroms)].copy()
    fr["x"] = fr["cb"].map(obs2x)
    fr["y"] = fr["cb"].map(obs2y)
    fr = fr.dropna(subset=["x", "y"])
    fr["x"] = fr["x"].astype(int)
    fr["y"] = fr["y"].astype(int)
    fr = fr[(fr.x >= 1) & (fr.x <= grid) & (fr.y >= 1) & (fr.y <= grid)]

    mid = (fr["start"].values + fr["end"].values) // 2
    coff = fr["chr"].map(off).values
    cmax = fr["chr"].map(lambda c: (sizes[c] + bin_size - 1) // bin_size - 1).values
    bidx = (coff + np.minimum(mid // bin_size, cmax)).astype(np.int64)
    sidx = ((fr["x"].values - 1) * grid + (fr["y"].values - 1)).astype(np.int64)

    n_spots = grid * grid
    counts = np.bincount(sidx * n_bins + bidx, minlength=n_spots * n_bins)
    counts = counts.reshape(n_spots, n_bins).astype(np.float32)

    xs = np.repeat(np.arange(1, grid + 1), grid)
    ys = np.tile(np.arange(1, grid + 1), grid)
    adata = _assemble(counts, xs, ys, bins, sample, bin_size, genome, grid)
    # drop empty spots (no fragments) up front — keeps downstream sparse-aware
    adata = adata[adata.obs["total_frags"] > 0].copy()
    return adata


def from_matrix(
    counts: np.ndarray,
    bins: pd.DataFrame,
    *,
    xs=None,
    ys=None,
    sample: str = "sample",
    bin_size: int = 1_000_000,
    genome: str = "hg38",
    grid: int = GRID,
) -> ad.AnnData:
    """Build an AnnData directly from an in-memory spot x bin count matrix.

    Used by :mod:`spatial_omics.datasets` and by anyone bringing their own matrix
    (e.g. the legacy ``*_counts.npz`` from create_count_matrix.py).
    """
    counts = np.asarray(counts, dtype=np.float32)
    if xs is None or ys is None:
        xs = np.repeat(np.arange(1, grid + 1), grid)
        ys = np.tile(np.arange(1, grid + 1), grid)
    return _assemble(counts, np.asarray(xs), np.asarray(ys), bins, sample, bin_size, genome, grid)


def _assemble(counts, xs, ys, bins, sample, bin_size, genome, grid) -> ad.AnnData:
    obs = pd.DataFrame({
        "x_id": np.asarray(xs, dtype=int),
        "y_id": np.asarray(ys, dtype=int),
        "total_frags": counts.sum(1).astype(int),
    })
    obs.index = [f"spot_{x}_{y}" for x, y in zip(obs.x_id, obs.y_id)]
    var = bins.copy()
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.obsm["spatial"] = np.column_stack([obs.x_id.values, obs.y_id.values]).astype(float)
    adata.uns["spatial_omics"] = {
        "sample": sample, "grid": int(grid), "bin_size": int(bin_size),
        "genome": genome, "version": "0.1.0",
    }
    return adata


def read_h5ad(path: str) -> ad.AnnData:
    """Read a spatial_omics object (thin wrapper around anndata.read_h5ad)."""
    return ad.read_h5ad(path)
