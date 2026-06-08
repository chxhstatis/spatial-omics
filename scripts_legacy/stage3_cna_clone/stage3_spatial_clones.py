#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage3_spatial_clones.py  (v3 — GC-corrected + paper tricks)

Stage-3 spatial clone analysis for spatial DNA-seq (DBiT 50x50, hg38, PDAC).
NO matched normal control -> in-data pseudo-normal (flattest cluster).

Pipeline: spot x 1Mb-bin matrix -> bin QC -> per-spot GC-quantile correction
-> library norm -> per-bin pseudo-normal -> per-spot center -> depth regress
-> spatial+genomic smooth -> PCA.

Incorporates two tricks from the slide-DNA-seq paper (Zhao et al. Nature 2022):
  * DOUBLE SMOOTHING of PC scores in BOTH PC-space and xy-space before KMeans
    (paper's pc_scores_smo_both) — denoises sparse spatial data without blurring
    real clone boundaries.
  * PERMUTATION Z-SCORE spatial significance per genomic region (paper Fig 2c):
    shuffle fragments->spots preserving per-spot totals, build an empirical null,
    z=(obs-mean)/std, subtract the normal-cluster mean -> signed -log10 p map.
"""

import argparse
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, to_rgb
from scipy.ndimage import uniform_filter, median_filter
from scipy.stats import norm as _norm
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import calinski_harabasz_score

HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895,
}
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
PDAC_GENES = {"KRAS": ("chr12", 25.21), "CDKN2A": ("chr9", 21.97),
              "TP53": ("chr17", 7.67), "SMAD4": ("chr18", 51.03),
              "MYC": ("chr8", 127.74)}
GRID = 50


def build_bins(bin_size):
    rows, offset, gi = [], {}, 0
    for c in CHROMS:
        n = (HG38[c] + bin_size - 1) // bin_size
        offset[c] = gi
        for b in range(n):
            rows.append((c, b * bin_size, min((b + 1) * bin_size, HG38[c])))
        gi += n
    return pd.DataFrame(rows, columns=["chr", "start", "end"]), offset, gi


def load_gc(gc_path, bins, bin_size):
    """align the per-bin GC table to the same bin order."""
    gc = pd.read_csv(gc_path, sep="\t")
    key = bins["chr"] + ":" + (bins["start"] // bin_size).astype(str)
    gmap = dict(zip(gc["chr"] + ":" + (gc["start"] // bin_size).astype(str), gc["gc"]))
    return key.map(gmap).fillna(0.0).values


def load_matrix(frag_path, matched_path, bin_size):
    m = pd.read_csv(matched_path, sep="\t",
                    usecols=["barcode_obs", "x_id", "y_id", "full_status"])
    m = m[m["full_status"] == "mapped"].dropna(subset=["x_id", "y_id"])
    m["x_id"] = m["x_id"].astype(int)
    m["y_id"] = m["y_id"].astype(int)
    obs2x = dict(zip(m.barcode_obs, m.x_id))
    obs2y = dict(zip(m.barcode_obs, m.y_id))

    bins, offset, n_bins = build_bins(bin_size)
    if frag_path.endswith(".gz"):
        proc = subprocess.Popen(["gzip", "-dc", frag_path], stdout=subprocess.PIPE)
        fr = pd.read_csv(proc.stdout, sep="\t", header=None,
                         names=["chr", "start", "end", "cb", "mapq"], on_bad_lines="skip")
        proc.stdout.close(); proc.wait()
    else:
        fr = pd.read_csv(frag_path, sep="\t", header=None,
                         names=["chr", "start", "end", "cb", "mapq"], on_bad_lines="skip")
    fr = fr[fr["chr"].isin(CHROMS)].copy()
    fr["x"] = fr["cb"].map(obs2x)
    fr["y"] = fr["cb"].map(obs2y)
    fr = fr.dropna(subset=["x", "y"])
    fr["x"] = fr["x"].astype(int); fr["y"] = fr["y"].astype(int)
    fr = fr[(fr.x >= 1) & (fr.x <= GRID) & (fr.y >= 1) & (fr.y <= GRID)]

    mid = (fr["start"].values + fr["end"].values) // 2
    coff = fr["chr"].map(offset).values
    cmax = fr["chr"].map(lambda c: (HG38[c] + bin_size - 1) // bin_size - 1).values
    local = np.minimum(mid // bin_size, cmax)
    bin_idx = (coff + local).astype(np.int64)
    spot_idx = ((fr["x"].values - 1) * GRID + (fr["y"].values - 1)).astype(np.int64)
    n_spots = GRID * GRID
    flat = spot_idx * n_bins + bin_idx
    counts = np.bincount(flat, minlength=n_spots * n_bins).reshape(n_spots, n_bins)
    return counts, bins, n_bins


def process(counts, gc, spatial_win, min_bin_frac=0.2, n_gc=20):
    """bin-QC -> spatial-smooth raw counts -> per-spot GC-quantile correction
       -> library norm -> per-bin reference."""
    tot = counts.sum(1).astype(float)
    pres = tot > 0
    bin_mean = counts[pres].mean(0)
    thr = min_bin_frac * np.median(bin_mean[bin_mean > 0])
    keep = bin_mean > thr

    C = counts[:, keep].astype(float)
    G = gc[keep]
    # densify before GC fit: spatially smooth raw counts on the grid
    C = spatial_smooth(C, pres, spatial_win)

    # GC-quantile groups (equal #bins per group by GC rank)
    ranks = np.argsort(np.argsort(G))
    grp = np.minimum((ranks * n_gc // len(G)), n_gc - 1)
    onehot = np.zeros((len(G), n_gc)); onehot[np.arange(len(G)), grp] = 1.0
    gcount = onehot.sum(0)                       # bins per GC group
    gsum = C @ onehot                            # (spots, n_gc) summed counts
    expected = np.divide(gsum, gcount[None, :], out=np.zeros_like(gsum), where=gcount[None, :] > 0)
    exp_bin = expected[:, grp]                   # per-spot expected count given GC
    rowmean = C.mean(1, keepdims=True)
    Cc = np.divide(C, exp_bin, out=np.zeros_like(C), where=exp_bin > 0) * rowmean

    # library-size normalize per spot
    t2 = Cc.sum(1)
    med = np.median(t2[pres]) if pres.any() else 1.0
    libn = np.zeros_like(Cc)
    libn[pres] = Cc[pres] / t2[pres, None] * med

    # per-bin pseudo-normal reference (cross-spot median)
    ref = np.median(libn[pres], axis=0); ref[ref <= 0] = np.nan
    rel = np.full_like(libn, np.nan)
    rel[pres] = libn[pres] / ref
    return rel, keep, tot, pres


def spatial_smooth(mat, present, win):
    n_bins = mat.shape[1]
    data = np.zeros((GRID, GRID, n_bins)); mask = np.zeros((GRID, GRID))
    filled = np.nan_to_num(mat, nan=0.0)
    for s in range(GRID * GRID):
        if present[s]:
            data[s // GRID, s % GRID] = filled[s]; mask[s // GRID, s % GRID] = 1.0
    num = uniform_filter(data * mask[:, :, None], size=(win, win, 1), mode="nearest")
    den = uniform_filter(mask, size=(win, win), mode="nearest")[:, :, None]
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0).reshape(GRID * GRID, n_bins)


def genome_smooth(L, bins, keep_bins, win=3):
    kb = bins[keep_bins].reset_index(drop=True)
    out = L.copy()
    for c in CHROMS:
        idx = np.where(kb["chr"].values == c)[0]
        if len(idx) >= win:
            out[:, idx] = median_filter(L[:, idx], size=(1, win), mode="nearest")
    return out


def choose_k(X, kmin, kmax, seed):
    scores, best = {}, (kmin, -1, None)
    for k in range(kmin, kmax + 1):
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
        s = calinski_harabasz_score(X, km.labels_)
        scores[k] = s
        if s > best[1]:
            best = (k, s, km.labels_)
    return best[0], best[2], scores


def double_smooth_pcs(pcs, xy_coords, k_pc=50, k_xy=50, k_join=10):
    """Paper trick (pc_scores_smo_both): smooth PC scores in PC-space, then average
    each spot's PC-smoothed value over its nearest xy neighbours. Denoises sparse
    spatial data for clustering without blurring real spatial clone boundaries."""
    n = len(pcs)
    kpc = min(k_pc, n); kxy = min(k_xy, n)
    _, idx_pc = NearestNeighbors(n_neighbors=kpc).fit(pcs).kneighbors(pcs)
    smo_pc = pcs[idx_pc].mean(axis=1)
    _, idx_xy = NearestNeighbors(n_neighbors=kxy).fit(xy_coords).kneighbors(xy_coords)
    kj = min(k_join, idx_xy.shape[1])
    return smo_pc[idx_xy[:, :kj]].mean(axis=1)          # joint PC+xy smoothing


def permutation_zscore_region(region_counts, spot_total, present, normal_mask,
                              xy_coords, n_perm, k_xy, seed):
    """Paper Fig 2c spatial significance: per-spot fragment count in a genomic
    region vs an empirical null built by randomly reassigning fragments to spots
    (multinomial, probability ~ each spot's share of total fragments), preserving
    the region's total. xy-smooth (sum over k neighbours), z-score, subtract the
    normal-cluster mean. Returns per-spot signed z (NaN where absent)."""
    idx = np.where(present)[0]
    obs = region_counts[idx].astype(float)
    tot = spot_total[idx].astype(float)
    M = int(round(obs.sum()))
    p = tot / tot.sum() if tot.sum() > 0 else np.ones(len(tot)) / len(tot)
    rng = np.random.default_rng(seed)
    perm = rng.multinomial(M, p, size=n_perm).astype(float)      # (n_perm, n_present)
    kk = min(k_xy, len(idx))
    _, nn = NearestNeighbors(n_neighbors=kk).fit(xy_coords).kneighbors(xy_coords)
    obs_s = obs[nn].sum(axis=1)
    perm_s = perm[:, nn].sum(axis=2)                             # (n_perm, n_present)
    mu = perm_s.mean(axis=0); sd = perm_s.std(axis=0) + 1e-9
    z = (obs_s - mu) / sd
    nm = normal_mask[idx]
    if nm.any():
        z = z - np.nanmean(z[nm])
    full = np.full(GRID * GRID, np.nan); full[idx] = z
    return full


def chrom_boundaries(bins, keep_bins):
    kb = bins[keep_bins].reset_index(drop=True)
    bnds, centers, order, pos = [], [], [], 0
    for c in CHROMS:
        n = int((kb["chr"].values == c).sum())
        if n == 0:
            continue
        bnds.append(pos); centers.append(pos + n / 2); order.append(c.replace("chr", "")); pos += n
    bnds.append(pos)
    return bnds, centers, order


def run_sample(sample, frag_path, matched_path, gc_path, outdir, bin_size,
               spatial_win, genome_win, kmin, kmax, n_pcs, seed,
               zscore_chroms=("chr8", "chr17", "chr18"), n_perm=100):
    os.makedirs(outdir, exist_ok=True)
    print("[%s] build %d-bp matrix" % (sample, bin_size))
    counts, bins, n_bins = load_matrix(frag_path, matched_path, bin_size)
    gc = load_gc(gc_path, bins, bin_size)

    rel, keep_bins, tot, present = process(counts, gc, spatial_win)
    kb = bins[keep_bins].reset_index(drop=True)
    kb.to_csv(os.path.join(outdir, "bin_info.tsv"), sep="\t", index=False)
    print("[%s] spots=%d bins kept=%d/%d" % (sample, int(present.sum()), int(keep_bins.sum()), n_bins))

    L = np.log2(np.clip(rel, 1e-3, None))
    L = np.nan_to_num(L, nan=0.0)
    L = L - np.median(L, axis=1, keepdims=True)        # per-spot center (remove ploidy/DC offset)

    # regress out the depth confound: residualize each bin against log10(coverage)
    cov = np.log10(tot + 1.0)
    covp = cov[present] - cov[present].mean()
    denom = float((covp ** 2).sum())
    if denom > 0:
        Lp = L[present]
        beta = (Lp * covp[:, None]).sum(axis=0) / denom     # per-bin slope vs depth
        L[present] = Lp - covp[:, None] * beta[None, :]

    L = genome_smooth(L, bins, keep_bins, win=genome_win)
    L = np.clip(L, -2.0, 2.0)

    Lp = L[present]
    pcs = PCA(n_components=min(n_pcs, Lp.shape[1], Lp.shape[0] - 1), random_state=seed).fit_transform(Lp)
    # paper trick: double-smooth PC scores (PC-space + xy-space) before clustering
    idxp = np.where(present)[0]
    xy_coords = np.column_stack([idxp // GRID, idxp % GRID]).astype(float)
    pcs_smo = double_smooth_pcs(pcs, xy_coords)
    k, lab_p, scores = choose_k(pcs_smo, kmin, kmax, seed)
    labels = np.full(GRID * GRID, -1); labels[present] = lab_p

    # depth-confound QC: correlation of PC1 with log10 total fragments
    r_depth = float(np.corrcoef(pcs[:, 0], np.log10(tot[present] + 1))[0, 1])

    flat = {c: np.median(np.abs(Lp[lab_p == c])) for c in range(k)}
    normal_c = min(flat, key=flat.get)

    clone_mean = np.vstack([L[labels == c].mean(axis=0) for c in range(k)])
    cna = clone_mean - clone_mean[normal_c][None, :]
    cna = cna - np.median(cna, axis=1, keepdims=True)   # median-center each profile (kill residual DC)

    aber = {c: np.mean(np.abs(cna[c])) for c in range(k)}
    tumor = sorted([c for c in range(k) if c != normal_c], key=lambda c: -aber[c])
    name = {normal_c: "Normal"}
    for i, c in enumerate(tumor, 1):
        name[c] = "Clone%d" % i

    print("[%s] k=%d sizes=%s normal=%s depth_confound(r_PC1_vs_cov)=%.2f" %
          (sample, k, np.bincount(lab_p).tolist(), name[normal_c], r_depth))

    xs = np.repeat(np.arange(1, GRID + 1), GRID); ys = np.tile(np.arange(1, GRID + 1), GRID)
    df = pd.DataFrame({"x_id": xs, "y_id": ys, "total_frags": tot.astype(int),
                       "present": present.astype(int), "cluster": labels})
    df["clone_label"] = df["cluster"].map(lambda c: name.get(c, "NA") if c >= 0 else "empty")
    df.to_csv(os.path.join(outdir, "spot_clone_labels.tsv"), sep="\t", index=False)

    cna_df = pd.concat([kb, pd.DataFrame(cna.T, columns=[name[c] for c in range(k)])], axis=1)
    cna_df.to_csv(os.path.join(outdir, "clone_cna_matrix.tsv"), sep="\t", index=False)

    with open(os.path.join(outdir, "summary.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        f.write("sample\t%s\nbin_size\t%d\nspots_present\t%d\nbins_kept\t%d\n" %
                (sample, bin_size, int(present.sum()), int(keep_bins.sum())))
        f.write("median_frags_per_spot\t%.1f\nchosen_k\t%d\nnormal_cluster\t%d\n" %
                (np.median(tot[present]), k, normal_c))
        f.write("depth_confound_r_PC1_vs_logcov\t%.3f\n" % r_depth)
        for c in range(k):
            f.write("cluster_%d_label\t%s\ncluster_%d_n_spots\t%d\ncluster_%d_mean_cov\t%.0f\n" %
                    (c, name[c], c, int((labels == c).sum()), c,
                     tot[(labels == c)].mean() if (labels == c).any() else 0))
    with open(os.path.join(outdir, "kmeans_model_scores.tsv"), "w") as f:
        f.write("k\tcalinski_harabasz\n")
        for kk, ss in scores.items():
            f.write("%d\t%.2f\n" % (kk, ss))

    # ---- paper Fig 2c: permutation z-score spatial significance per region ----
    normal_mask = (labels == normal_c)
    zmaps = {}
    for chrom in zscore_chroms:
        region = (bins["chr"].values == chrom)
        if not region.any():
            continue
        rc = counts[:, region].sum(axis=1)
        z = permutation_zscore_region(rc, tot, present, normal_mask, xy_coords,
                                      n_perm, k_xy=10, seed=seed)
        zmaps[chrom] = z
    if zmaps:
        zdf = pd.DataFrame({"x_id": xs, "y_id": ys})
        for c, z in zmaps.items():
            zdf["z_" + c] = z
        zdf.to_csv(os.path.join(outdir, "region_zscores.tsv"), sep="\t", index=False)
        _plot_zmaps(sample, outdir, zmaps, present)

    _plot_all(sample, outdir, labels, name, k, normal_c, cna, kb, bins, keep_bins, tot, present)
    print("[%s] DONE -> %s" % (sample, outdir))


def _plot_zmaps(sample, outdir, zmaps, present):
    """signed -log10 p spatial maps (paper Fig 2c style: red=gain, blue=loss)."""
    chroms = list(zmaps.keys()); ncol = len(chroms)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 4), squeeze=False)
    for ax, c in zip(axes[0], chroms):
        z = zmaps[c]
        # signed -log10 two-sided p from z
        p = -np.log10(np.clip(2 * _norm.cdf(-np.abs(z)), 1e-300, 1.0)) * np.sign(z)
        grid = np.full((GRID, GRID), np.nan)
        for s in range(GRID * GRID):
            if present[s] and not np.isnan(p[s]):
                grid[s // GRID, s % GRID] = p[s]
        im = ax.imshow(grid, origin="lower", cmap="RdBu_r",
                       norm=TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5), interpolation="nearest")
        ax.set_title("%s  signed -log10 p" % c, fontsize=10)
        ax.set_xlabel("Y"); ax.set_ylabel("X")
        plt.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("%s: spatial significance of regional coverage (permutation z, paper Fig 2c style)" % sample,
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "region_significance_maps.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "region_significance_maps.pdf")); plt.close()


def _plot_all(sample, outdir, labels, name, k, normal_c, cna, kb, bins, keep_bins, tot, present):
    grid_lab = labels.reshape(GRID, GRID)
    grid_cov = np.where(present, tot, np.nan).reshape(GRID, GRID)
    clone_colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf"]
    colmap = {normal_c: "#9e9e9e"}; ci = 0
    for c in range(k):
        if c != normal_c:
            colmap[c] = clone_colors[ci % len(clone_colors)]; ci += 1
    rgb = np.ones((GRID, GRID, 3))
    for x in range(GRID):
        for y in range(GRID):
            c = grid_lab[x, y]
            rgb[x, y] = to_rgb(colmap[c]) if c >= 0 else (1, 1, 1)
    plt.figure(figsize=(7, 6.5))
    plt.imshow(rgb, origin="lower", interpolation="nearest")
    handles = [plt.Line2D([0], [0], marker="s", ls="", markersize=10, markerfacecolor=colmap[c],
               label="%s (n=%d)" % (name[c], int((labels == c).sum()))) for c in range(k)]
    plt.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    plt.xlabel("Y"); plt.ylabel("X"); plt.title("%s: spatial clone map (50x50, 50um)" % sample)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "spatial_clone_map.png"), dpi=220)
    plt.savefig(os.path.join(outdir, "spatial_clone_map.pdf")); plt.close()

    plt.figure(figsize=(7, 6))
    im = plt.imshow(np.log10(grid_cov), origin="lower", cmap="viridis", interpolation="nearest")
    plt.colorbar(im, label="log10(fragments per spot)")
    plt.xlabel("Y"); plt.ylabel("X"); plt.title("%s: fragments per spot" % sample)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "spot_coverage_map.png"), dpi=200); plt.close()

    bnds, centers, order = chrom_boundaries(bins, keep_bins)
    plt.figure(figsize=(14, 0.7 * k + 2))
    plt.imshow(cna, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), interpolation="nearest")
    plt.colorbar(label="log2 ratio vs in-sample normal (median-centered)")
    plt.yticks(range(k), [name[c] for c in range(k)])
    for b in bnds:
        plt.axvline(b - 0.5, color="k", lw=0.4, alpha=0.5)
    plt.xticks(centers, order, fontsize=7)
    plt.xlabel("Genomic position (chr)"); plt.title("%s: clone copy-number (GC-corrected)" % sample)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "clone_cna_heatmap.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "clone_cna_heatmap.pdf")); plt.close()

    tumor_cs = [c for c in range(k) if c != normal_c]
    if tumor_cs:
        kbr = kb.reset_index(drop=True)
        gene_x = {}
        for g, (gc_, gmb) in PDAC_GENES.items():
            sel = kbr.index[(kbr["chr"] == gc_) & (kbr["start"] <= gmb * 1e6) & (kbr["end"] > gmb * 1e6)]
            if len(sel):
                gene_x[g] = sel[0]
        fig, axes = plt.subplots(len(tumor_cs), 1, figsize=(14, 2.2 * len(tumor_cs)), squeeze=False)
        for ax, c in zip(axes[:, 0], tumor_cs):
            ax.plot(cna[c], lw=0.5, color="#333")
            ax.axhline(0, color="grey", lw=0.5); ax.set_ylim(-1.2, 1.2)
            for b in bnds:
                ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.4)
            for g, gx in gene_x.items():
                ax.axvline(gx, color="darkorange", lw=0.6, ls="--", alpha=0.8)
                ax.text(gx, 1.05, g, fontsize=7, rotation=90, va="bottom", color="darkorange")
            ax.set_xticks(centers); ax.set_xticklabels(order, fontsize=7); ax.set_ylabel("%s\nlog2" % name[c])
        axes[0, 0].set_title("%s: per-clone genome-wide CNA (GC-corrected, vs in-sample normal)" % sample)
        plt.tight_layout(); plt.savefig(os.path.join(outdir, "clone_cna_profiles.png"), dpi=200)
        plt.savefig(os.path.join(outdir, "clone_cna_profiles.pdf")); plt.close()


def resolve_sample_dir(stage2_dir, sample):
    """accept either <stage2_dir>/<sample>/ or HPC-style <stage2_dir>/<sample>_xy_map/."""
    for cand in (sample, sample + "_xy_map", sample + "_out"):
        d = os.path.join(stage2_dir, cand)
        if os.path.isfile(os.path.join(d, "top50x50_fragments_with_cb.tsv.gz")):
            return d
    return os.path.join(stage2_dir, sample)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage2_dir", required=True)
    ap.add_argument("--samples", nargs="+", default=["520_520", "525_525"])
    ap.add_argument("--gc_file", required=True)
    ap.add_argument("--outroot", required=True)
    ap.add_argument("--bin_size", type=int, default=1_000_000)
    ap.add_argument("--spatial_win", type=int, default=5)
    ap.add_argument("--genome_win", type=int, default=3)
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--n_pcs", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--zscore_chroms", default="chr8,chr17,chr18",
                    help="comma-separated chromosomes for permutation z-score spatial significance maps")
    ap.add_argument("--n_perm", type=int, default=100)
    args = ap.parse_args()
    zchr = tuple(c.strip() for c in args.zscore_chroms.split(",") if c.strip())
    for s in args.samples:
        sd = resolve_sample_dir(args.stage2_dir, s)
        frag = os.path.join(sd, "top50x50_fragments_with_cb.tsv.gz")
        matched = os.path.join(sd, "matched_spot_barcodes.tsv")
        if not os.path.isfile(frag):
            print("[%s] SKIP missing %s" % (s, frag)); continue
        try:
            run_sample(s, frag, matched, args.gc_file, os.path.join(args.outroot, s),
                       args.bin_size, args.spatial_win, args.genome_win,
                       args.kmin, args.kmax, args.n_pcs, args.seed,
                       zscore_chroms=zchr, n_perm=args.n_perm)
        except Exception as e:
            import traceback
            print("[%s] FAILED: %s" % (s, e)); traceback.print_exc()
    print("[ALL DONE]")


if __name__ == "__main__":
    main()
