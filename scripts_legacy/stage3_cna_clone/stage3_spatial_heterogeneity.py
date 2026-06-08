#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage3_spatial_heterogeneity.py

Per-sample intratumour spatial CNV heterogeneity, built on the per-spot CNV matrix
(stage3_per_spot_cnv.py output). The honest question first: IS there real spatial
structure, or is the spot-to-spot variation just (a) tissue-coverage density or (b) the
spatial smoothing kernel? We answer it with Moran's I, but against the right controls.

!!! TWO CONFOUNDS THIS TOOL CONTROLS FOR (do not skip) !!!
1. SMOOTHING: per-spot CNV is built with a spatial window (--spatial_win, default 7), so
   neighbouring spots SHARE input data and are similar BY CONSTRUCTION -> Moran's I on the
   smoothed matrix is inflated (0.7-0.9) and is NOT evidence of biology. Pass a win=1
   (unsmoothed) matrix via --raw_perspot_root for the honest burden autocorrelation.
2. COVERAGE: tissue density (fragments/spot) is itself spatially autocorrelated; CNV burden
   partly tracks coverage. We compute Moran's I of coverage as the baseline. Real CNA
   heterogeneity must EXCEED the coverage autocorrelation, measured unsmoothed.
The verdict 'cna_beyond_coverage' is True only if raw(win=1) burden Moran > coverage Moran.

Per sample:
  1. Moran's I of per-spot CNV burden (8-neighbour weights) + permutation p.
  2. PCA of the spot x bin CNV matrix; Moran's I of PC1/PC2 (spatial structure in the
     copy-number profiles themselves, not just total burden).
  3. KMeans on top PCs (k by silhouette) -> spatial CNV subregions mapped on the 50x50 grid.
  4. Per-subregion mean CNV profile (genome-wide).

Outputs (per sample):
  heterogeneity_summary.tsv   moran_burden/p, moran_pc1/p, moran_pc2/p, k, var_explained
  spatial_subregion_map.png   50x50 grid coloured by CNV subregion
  cnv_burden_map.png          50x50 grid of per-spot CNV burden
  subregion_cnv_profiles.png  mean copy number per subregion across the genome
  spot_subregion_labels.tsv   x_id,y_id,total_frags,subregion

A cohort summary (Moran's I per sample) is written to the output root.
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

GRID = 50
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]


def neighbours8():
    """8-neighbour adjacency over the 50x50 grid -> list of (i, j) index pairs (i<j)."""
    idx = {}
    for x in range(GRID):
        for y in range(GRID):
            idx[(x, y)] = x * GRID + y
    pairs = []
    for x in range(GRID):
        for y in range(GRID):
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < GRID and 0 <= ny < GRID:
                        a, b = idx[(x, y)], idx[(nx, ny)]
                        if a < b:
                            pairs.append((a, b))
    return np.array(pairs)


PAIRS = neighbours8()


def morans_I(values, present, n_perm=999, seed=0):
    """Moran's I on the 50x50 grid using 8-neighbour weights; values NaN where absent.
    Returns (I, p_perm). Only edges between two present spots count.

    Significance via a TRUE random permutation null: the observed spot values are
    randomly reassigned to spot positions (label permutation), which destroys spatial
    arrangement while preserving the value distribution. (An earlier version used a
    cyclic shift, which PRESERVES local autocorrelation and is therefore an invalid
    null — fixed here. The cna_beyond_coverage verdict never used this p-value; it is a
    direct comparison of two observed I values, so that verdict is unaffected.)"""
    v = values.copy().astype(float)
    ok = present & np.isfinite(v)
    edge = PAIRS[ok[PAIRS[:, 0]] & ok[PAIRS[:, 1]]]
    if len(edge) < 10 or ok.sum() < 20:
        return np.nan, np.nan
    z = np.zeros_like(v); m = v[ok].mean(); z[ok] = v[ok] - m
    denom = (z[ok] ** 2).sum()

    def I_of(zz):
        num = (zz[edge[:, 0]] * zz[edge[:, 1]]).sum()  # each undirected edge once
        return (ok.sum() / len(edge)) * (num / denom)

    I = I_of(z)
    rng = np.random.default_rng(seed)             # reproducible random permutations
    order = ok.nonzero()[0]
    vals = v[ok].copy()
    null = np.empty(n_perm)
    for p in range(n_perm):
        perm = rng.permutation(vals)              # genuine label permutation
        zz = np.zeros_like(v); zz[order] = perm - perm.mean()
        null[p] = I_of(zz)
    p_val = (1 + int((null >= I).sum())) / (1 + n_perm)
    return float(I), float(p_val)


def _vec(df, col):
    v = np.full(GRID * GRID, np.nan)
    v[((df.x_id.values - 1) * GRID + (df.y_id.values - 1)).astype(int)] = df[col].values.astype(float)
    return v


def run_sample(perspot_dir, sample, outdir, raw_perspot_dir=None):
    os.makedirs(outdir, exist_ok=True)
    M = pd.read_csv(os.path.join(perspot_dir, "per_spot_cnv_matrix.tsv.gz"), sep="\t")
    bincols = [c for c in M.columns if ":" in c]
    chrom_of = np.array([c.split(":")[0] for c in bincols])
    X = M[bincols].values.astype(float)
    sidx = ((M.x_id.values - 1) * GRID + (M.y_id.values - 1)).astype(int)
    present = np.zeros(GRID * GRID, bool); present[sidx] = True

    bdf = pd.read_csv(os.path.join(perspot_dir, "per_spot_cnv_burden.tsv"), sep="\t")
    bvec = _vec(bdf, "cnv_burden")

    # ---- controls: coverage autocorrelation, and raw(win=1) burden autocorrelation ----
    cov_vec = _vec(bdf, "total_frags")
    mi_cov, _ = morans_I(cov_vec, present)
    mi_b_raw = np.nan
    if raw_perspot_dir and os.path.isfile(os.path.join(raw_perspot_dir, "per_spot_cnv_burden.tsv")):
        rdf = pd.read_csv(os.path.join(raw_perspot_dir, "per_spot_cnv_burden.tsv"), sep="\t")
        rvec = _vec(rdf, "cnv_burden"); rpres = present & np.isfinite(rvec)
        if np.nanstd(rvec[rpres]) > 0:
            mi_b_raw, _ = morans_I(rvec, present)
    cna_beyond_cov = bool(np.isfinite(mi_b_raw) and np.isfinite(mi_cov) and mi_b_raw > mi_cov)

    # ---- 1. Moran's I of (smoothed) burden — INFLATED by smoothing, reported with caveat ----
    mi_b, p_b = morans_I(bvec, present)

    # ---- 2. PCA + Moran's I of PCs ----
    Xc = np.nan_to_num(X - np.nanmean(X, 0))
    npc = min(10, Xc.shape[0] - 1, Xc.shape[1])
    pca = PCA(n_components=npc, random_state=0).fit(Xc)
    pcs = pca.transform(Xc)  # rows aligned to M
    pc_grid = np.full((GRID * GRID, npc), np.nan); pc_grid[sidx] = pcs
    mi_pc, p_pc = [], []
    for k in range(min(3, npc)):
        I, p = morans_I(pc_grid[:, k], present)
        mi_pc.append(I); p_pc.append(p)

    # ---- 3. KMeans subregions (k by silhouette on top PCs) ----
    top = pcs[:, :min(5, npc)]
    best_k, best_s, best_lab = 2, -1, None
    for k in range(2, 7):
        lab = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(top)
        try:
            sc = silhouette_score(top, lab)
        except Exception:
            sc = -1
        if sc > best_s:
            best_k, best_s, best_lab = k, sc, lab
    lab_grid = np.full(GRID * GRID, -1); lab_grid[sidx] = best_lab

    # ---- outputs ----
    pd.DataFrame({"x_id": M.x_id, "y_id": M.y_id, "total_frags": M.total_frags,
                  "subregion": best_lab}).to_csv(os.path.join(outdir, "spot_subregion_labels.tsv"), sep="\t", index=False)
    summ = dict(sample=sample, n_spots=int(present.sum()), n_bins=len(bincols),
                moran_coverage=round(mi_cov, 3),
                moran_burden_raw_win1=round(mi_b_raw, 3) if np.isfinite(mi_b_raw) else np.nan,
                moran_burden_smoothed=round(mi_b, 3),
                cna_beyond_coverage=cna_beyond_cov,
                moran_pc1_smoothed=round(mi_pc[0], 3),
                k=best_k, silhouette=round(best_s, 3),
                pc1_var=round(float(pca.explained_variance_ratio_[0]), 3))
    pd.DataFrame([summ]).to_csv(os.path.join(outdir, "heterogeneity_summary.tsv"), sep="\t", index=False)

    # burden map
    g = np.full((GRID, GRID), np.nan)
    for s in np.where(present)[0]:
        g[s // GRID, s % GRID] = bvec[s]
    plt.figure(figsize=(6.5, 5.5))
    im = plt.imshow(g, origin="lower", cmap="magma", interpolation="nearest")
    plt.colorbar(im, label="CNV burden"); plt.xlabel("Y"); plt.ylabel("X")
    plt.title("%s: CNV burden | coverage I=%.2f, raw I=%s, smoothed I=%.2f (smoothed inflated)"
              % (sample, mi_cov, ("%.2f" % mi_b_raw if np.isfinite(mi_b_raw) else "flat"), mi_b))
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "cnv_burden_map.png"), dpi=200); plt.close()

    # subregion map
    gl = np.full((GRID, GRID), np.nan)
    for s in np.where(present)[0]:
        gl[s // GRID, s % GRID] = lab_grid[s]
    plt.figure(figsize=(6.5, 5.5))
    im = plt.imshow(gl, origin="lower", cmap="tab10", interpolation="nearest", vmin=0, vmax=9)
    plt.colorbar(im, label="CNV subregion"); plt.xlabel("Y"); plt.ylabel("X")
    plt.title("%s: CNV subregions (k=%d, silhouette=%.2f)" % (sample, best_k, best_s))
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "spatial_subregion_map.png"), dpi=200); plt.close()

    # per-subregion mean profile
    order_idx = np.argsort([CHROMS.index(c) if c in CHROMS else 99 for c in chrom_of])
    bnds, centers, labs, pos = [], [], [], 0
    co = chrom_of[order_idx]
    for c in CHROMS:
        m = int((co == c).sum())
        if m: bnds.append(pos); centers.append(pos + m / 2); labs.append(c.replace("chr", "")); pos += m
    fig, ax = plt.subplots(figsize=(14, 4))
    for cl in range(best_k):
        prof = np.nanmean(X[best_lab == cl], 0)[order_idx]
        ax.plot(prof, lw=0.8, label="subregion %d (n=%d)" % (cl, int((best_lab == cl).sum())))
    ax.axhline(2, color="grey", lw=0.5)
    for b in bnds: ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.3)
    ax.set_xticks(centers); ax.set_xticklabels(labs, fontsize=7); ax.set_ylim(1, 3)
    ax.set_ylabel("mean copy number"); ax.legend(fontsize=7, ncol=2)
    ax.set_title("%s: mean CNV profile per spatial subregion" % sample)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "subregion_cnv_profiles.png"), dpi=200); plt.close()

    verdict = "CNA>coverage (real)" if cna_beyond_cov else "NOT beyond coverage/smoothing"
    print("[%s] Moran: coverage=%.2f  burden_raw(win1)=%s  burden_smoothed=%.2f  -> %s"
          % (sample, mi_cov, ("%.2f" % mi_b_raw if np.isfinite(mi_b_raw) else "flat"), mi_b, verdict))
    return summ


def main():
    ap = argparse.ArgumentParser(description="Per-sample intratumour spatial CNV heterogeneity.")
    ap.add_argument("--perspot_root", required=True, help="stage3_out_perspot dir (holds <sample>/per_spot_cnv_matrix.tsv.gz)")
    ap.add_argument("--raw_perspot_root", default=None,
                    help="win=1 (unsmoothed) per-spot output root -> honest burden autocorrelation & verdict")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    rows = []
    for s in a.samples:
        d = os.path.join(a.perspot_root, s)
        if not os.path.isfile(os.path.join(d, "per_spot_cnv_matrix.tsv.gz")):
            print("[%s] SKIP (no per_spot_cnv_matrix)" % s); continue
        rawd = os.path.join(a.raw_perspot_root, s) if a.raw_perspot_root else None
        try:
            rows.append(run_sample(d, s, os.path.join(a.out, s), raw_perspot_dir=rawd))
        except Exception as ex:
            print("[%s] FAILED: %s" % (s, ex))
    if rows:
        summ = pd.DataFrame(rows)
        summ.to_csv(os.path.join(a.out, "cohort_heterogeneity_summary.tsv"), sep="\t", index=False)
        # bar chart: coverage vs raw-burden vs smoothed-burden (shows smoothing inflation + coverage baseline)
        fig, ax = plt.subplots(figsize=(8, 4)); x = np.arange(len(summ)); w = 0.27
        ax.bar(x - w, summ.moran_coverage, w, label="coverage (tissue density)", color="#7f8c8d")
        ax.bar(x, summ.moran_burden_raw_win1, w, label="CNV burden raw (win1)", color="#8e44ad")
        ax.bar(x + w, summ.moran_burden_smoothed, w, label="CNV burden smoothed (inflated)", color="#d2b4de")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(x); ax.set_xticklabels(summ["sample"], rotation=30, ha="right")
        ax.set_ylabel("Moran's I"); ax.legend(fontsize=8)
        ax.set_title("Spatial structure: real CNA heterogeneity (raw) must exceed coverage baseline")
        plt.tight_layout(); plt.savefig(os.path.join(a.out, "cohort_morans_I.png"), dpi=200); plt.close()
        print("\n[cohort] summary:\n", summ.to_string(index=False))


if __name__ == "__main__":
    main()
