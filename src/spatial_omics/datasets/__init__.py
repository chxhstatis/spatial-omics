"""Synthetic & public datasets for tutorials, tests and benchmarking.

``simulate()`` generates a fully synthetic spatial DNA-seq experiment with a known
ground-truth clone layout, so every tutorial and test runs with no external data
(patient data is never shared). It models the things that make real sparse spatial
DNA-seq hard: ~0.x genome coverage, GC/mappability bias, and spatial clone domains.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from .._bins import build_bins
from .._constants import chrom_order
from ..io import from_matrix

__all__ = ["simulate", "simulate_rna"]


def simulate(
    *,
    grid: int = 50,
    bin_size: int = 1_000_000,
    genome: str = "hg38",
    coverage: float = 0.4,
    n_clones: int = 2,
    seed: int = 0,
    sample: str = "sim",
):
    """Simulate a DBiT spatial DNA-seq sample with known clone structure.

    A central tumour region is split into ``n_clones`` spatial domains, each with a
    distinct copy-number profile (gains/losses on a few chromosomes); the surrounding
    tissue is diploid 'Normal'. Fragment counts are Poisson-drawn from
    ``copy_number x gc_bias x mappability``, scaled so the mean per-spot depth matches
    ``coverage`` (fraction of a 1x genome). Ground truth -> ``obs['true_clone']`` and
    ``uns['spatial_omics_sim']['true_profiles']``.

    Returns
    -------
    AnnData with ``layers['counts']``, ``var[gc/mappability/blacklist]``, ground truth.
    """
    rng = np.random.default_rng(seed)
    bins = build_bins(bin_size, genome)
    n_bins = len(bins)
    chroms = chrom_order(genome)

    # --- per-bin bias tracks (plausible-looking, bounded) ---
    gc = np.clip(rng.normal(0.45, 0.06, n_bins), 0.25, 0.65)
    mappability = np.clip(rng.beta(8, 1.2, n_bins), 0.3, 1.0)
    blacklist = np.where(rng.random(n_bins) < 0.03, rng.uniform(0.3, 1.0, n_bins), 0.0)
    bins = bins.assign(gc=gc, mappability=mappability, blacklist=blacklist)

    # --- ground-truth clone layout on the grid ---
    xs = np.repeat(np.arange(1, grid + 1), grid)
    ys = np.tile(np.arange(1, grid + 1), grid)
    cx, cy = (grid + 1) / 2, (grid + 1) / 2
    radius = grid * 0.32
    in_tumor = ((xs - cx) ** 2 + (ys - cy) ** 2) < radius ** 2
    angle = np.arctan2(ys - cy, xs - cx)  # split tumour into angular wedges = clones
    wedge = ((angle + np.pi) / (2 * np.pi) * n_clones).astype(int).clip(0, n_clones - 1)
    true_clone = np.where(in_tumor, np.array([f"Clone{w + 1}" for w in wedge]), "Normal")

    # --- copy-number profile per clone ---
    def _profile(events):
        cn = np.full(n_bins, 2.0)
        for chrom, lo, hi, val in events:
            m = (bins["chr"].values == chrom) & (bins["start"].values >= lo * 1e6) \
                & (bins["start"].values < hi * 1e6)
            cn[m] = val
        return cn

    # PDAC-flavoured events: chr9 (CDKN2A) & chr17 (TP53) & chr18 (SMAD4) losses, chr8 (MYC) gain
    clone_events = {
        "Clone1": [("chr9", 0, 40, 1), ("chr17", 0, 25, 1), ("chr8", 100, 146, 3)],
        "Clone2": [("chr18", 0, 60, 1), ("chr8", 100, 146, 4), ("chr12", 0, 40, 3)],
        "Clone3": [("chr9", 0, 40, 1), ("chr18", 0, 60, 1), ("chr20", 0, 64, 3)],
    }
    profiles = {"Normal": np.full(n_bins, 2.0)}
    for c in range(1, n_clones + 1):
        name = f"Clone{c}"
        profiles[name] = _profile(clone_events.get(name, []))

    # --- expected per-(spot,bin) rate ---
    gc_bias = np.exp(-((gc - 0.45) ** 2) / (2 * 0.08 ** 2))  # peaked GC efficiency
    bin_eff = gc_bias * mappability * (blacklist < 0.3)
    bin_eff = bin_eff / bin_eff.sum()

    # target mean fragments per spot for the requested coverage
    # (1x genome over n_bins ~ n_bins reads/bin baseline; scale by coverage)
    mean_frags = max(coverage * n_bins, 50.0)
    depth = rng.gamma(shape=4.0, scale=mean_frags / 4.0, size=grid * grid)  # spot depth variation

    counts = np.zeros((grid * grid, n_bins), dtype=np.float32)
    for name in profiles:
        sel = np.where(true_clone == name)[0]
        if len(sel) == 0:
            continue
        rate = (profiles[name] / 2.0) * bin_eff  # copy number scales expected coverage
        rate = rate / rate.sum()
        lam = depth[sel][:, None] * rate[None, :]
        counts[sel] = rng.poisson(lam).astype(np.float32)

    adata = from_matrix(counts, bins, xs=xs, ys=ys, sample=sample,
                        bin_size=bin_size, genome=genome, grid=grid)
    adata.var["gc"] = gc
    adata.var["mappability"] = mappability
    adata.var["blacklist"] = blacklist
    adata.obs["true_clone"] = pd.Categorical(true_clone)
    adata.uns["spatial_omics_sim"] = {
        "coverage": coverage, "n_clones": n_clones, "seed": seed,
        "true_profiles": pd.DataFrame(profiles, index=adata.var_names),
    }
    return adata[adata.obs["total_frags"] > 0].copy()


def simulate_rna(dna, *, n_genes=200, frac_genetic=0.3, frac_microenv=0.3, seed=1):
    """Simulate matched spatial RNA on the same grid as a ``simulate()`` DNA sample.

    Cell types are laid out so microenvironment is spatially structured: ``Tumor`` in
    the tumour disc, ``Immune`` in a margin ring around it, ``Normal`` elsewhere. Then:

    * ``frac_genetic`` of genes have a per-clone mean shift (intrinsic / genetic driver);
    * ``frac_microenv`` are driven by local tumour density (extrinsic / microenvironment);
    * the rest are pure noise.

    Ground-truth driver per gene -> ``var['true_driver']``. Pairs with
    :func:`spatial_omics.tl.variance_decomposition`, which should recover the split.
    """
    rng = np.random.default_rng(seed)
    xs = dna.obs["x_id"].values.astype(float)
    ys = dna.obs["y_id"].values.astype(float)
    n = len(xs)
    grid = int(dna.uns.get("spatial_omics", {}).get("grid", 50))
    cx = cy = (grid + 1) / 2
    r = np.hypot(xs - cx, ys - cy)
    r_tumor = grid * 0.32
    cell_type = np.where(r < r_tumor, "Tumor",
                         np.where(r < r_tumor + 4, "Immune", "Normal"))

    xy = np.column_stack([xs, ys])
    kk = min(50, n)
    nn = NearestNeighbors(n_neighbors=kk).fit(xy).kneighbors(xy, return_distance=False)
    tumor_density = (cell_type[nn] == "Tumor").mean(axis=1)

    # clone identity for the genetic driver (use called or true clone if present, else cell type)
    clone_src = dna.obs["clone"] if "clone" in dna.obs else dna.obs.get("true_clone")
    clones = clone_src.astype(str).values if clone_src is not None else cell_type
    uniq = pd.unique(clones)
    clone_code = {c: i for i, c in enumerate(uniq)}
    clone_idx = np.array([clone_code[c] for c in clones])

    n_gen = int(n_genes * frac_genetic)
    n_mic = int(n_genes * frac_microenv)
    drivers = (["genetic"] * n_gen + ["microenv"] * n_mic
               + ["noise"] * (n_genes - n_gen - n_mic))
    rng.shuffle(drivers)

    X = np.zeros((n, n_genes), dtype=np.float32)
    clone_effect = rng.normal(0, 2.0, (len(uniq), n_genes))  # per-clone mean per gene
    td = (tumor_density - tumor_density.mean()) / (tumor_density.std() + 1e-9)
    for g, drv in enumerate(drivers):
        base = rng.normal(5.0, 0.5)
        noise = rng.normal(0, 1.0, n)
        if drv == "genetic":
            X[:, g] = base + clone_effect[clone_idx, g] + noise
        elif drv == "microenv":
            X[:, g] = base + rng.uniform(2, 4) * td + noise
        else:
            X[:, g] = base + noise

    obs = pd.DataFrame({"x_id": dna.obs["x_id"].values, "y_id": dna.obs["y_id"].values,
                        "cell_type": pd.Categorical(cell_type)})
    obs.index = dna.obs_names
    var = pd.DataFrame({"true_driver": drivers}, index=[f"gene_{i}" for i in range(n_genes)])
    rna = ad.AnnData(X=X, obs=obs, var=var)
    rna.obsm["spatial"] = xy
    rna.uns["spatial_omics_rna_sim"] = {"frac_genetic": frac_genetic, "frac_microenv": frac_microenv}
    return rna
