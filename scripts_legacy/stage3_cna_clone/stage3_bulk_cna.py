#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage3_bulk_cna.py  (deliverable D)

For the signal-weak case (low PDAC purity, no robust spatial subclones), produce:
  1) a clean BULK tumor CNA profile (GC + mappability corrected, blacklist removed,
     normalized to modal=diploid) with PDAC genes annotated;
  2) a SPATIAL CNA-BURDEN / tumor-likeness map = per-spot similarity to the bulk
     tumor CNA profile, plotted on the 50x50 grid (a tumor-content proxy).

Inputs (per sample) from Stage-2: top50x50_fragments_with_cb.tsv.gz + matched_spot_barcodes.tsv
Reference: bin_gc_1mb.tsv, k100.umap.bed.gz (mappability), hg38-blacklist.v2.bed.gz
"""
import argparse
import gzip
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter, median_filter

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
PDAC_GENES = {"KRAS": ("chr12", 25.21), "CDKN2A": ("chr9", 21.97), "TP53": ("chr17", 7.67),
              "SMAD4": ("chr18", 51.03), "MYC": ("chr8", 127.74), "GATA6": ("chr18", 22.17)}
GRID = 50
BIN = 1_000_000


def load_normal_mask(path):
    """read a normal/reference spot list (cols x_id,y_id) -> boolean mask over GRID*GRID, or None."""
    if not path or not os.path.isfile(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if "x_id" not in df or "y_id" not in df:
        return None
    sidx = ((df.x_id.astype(int) - 1) * GRID + (df.y_id.astype(int) - 1)).values
    sidx = sidx[(sidx >= 0) & (sidx < GRID * GRID)]
    m = np.zeros(GRID * GRID, bool); m[sidx] = True
    return m


def build_bins():
    rows, off, gi = [], {}, 0
    for c in CHROMS:
        n = (HG38[c] + BIN - 1) // BIN
        off[c] = gi
        for b in range(n):
            rows.append((c, b * BIN, min((b + 1) * BIN, HG38[c])))
        gi += n
    return pd.DataFrame(rows, columns=["chr", "start", "end"]), off, gi


def bed_frac_per_bin(path, off, n_bins, skip_track):
    cov = np.zeros(n_bins)
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with op as f:
        for line in f:
            if skip_track and line.startswith("track"):
                continue
            p = line.split("\t")
            if len(p) < 3 or p[0] not in CHROMS:
                continue
            c, s, e = p[0], int(p[1]), int(p[2])
            b0, b1 = s // BIN, (e - 1) // BIN
            for b in range(b0, b1 + 1):
                ss, ee = max(s, b * BIN), min(e, (b + 1) * BIN)
                cov[off[c] + b] += (ee - ss)
    return cov


def load_counts(frag_path, matched_path, off, n_bins):
    m = pd.read_csv(matched_path, sep="\t", usecols=["barcode_obs", "x_id", "y_id", "full_status"])
    m = m[m.full_status == "mapped"].dropna(subset=["x_id", "y_id"])
    obs2x = dict(zip(m.barcode_obs, m.x_id.astype(int)))
    obs2y = dict(zip(m.barcode_obs, m.y_id.astype(int)))
    proc = subprocess.Popen(["gzip", "-dc", frag_path], stdout=subprocess.PIPE)
    fr = pd.read_csv(proc.stdout, sep="\t", header=None,
                     names=["chr", "start", "end", "cb", "mapq"], on_bad_lines="skip")
    proc.stdout.close(); proc.wait()
    fr = fr[fr.chr.isin(CHROMS)]
    fr["x"] = fr.cb.map(obs2x); fr["y"] = fr.cb.map(obs2y)
    fr = fr.dropna(subset=["x", "y"])
    fr["x"] = fr.x.astype(int); fr["y"] = fr.y.astype(int)
    fr = fr[(fr.x >= 1) & (fr.x <= GRID) & (fr.y >= 1) & (fr.y <= GRID)]
    mid = (fr.start.values + fr.end.values) // 2
    coff = fr.chr.map(off).values
    cmax = fr.chr.map(lambda c: (HG38[c] + BIN - 1) // BIN - 1).values
    bidx = (coff + np.minimum(mid // BIN, cmax)).astype(np.int64)
    sidx = ((fr.x.values - 1) * GRID + (fr.y.values - 1)).astype(np.int64)
    counts = np.bincount(sidx * n_bins + bidx, minlength=GRID * GRID * n_bins).reshape(GRID * GRID, n_bins)
    return counts


def quantile_correct(vec, covar, keep, n_q=30):
    """flatten vec vs covar by quantile-median; returns corrected vec (only on `keep`)."""
    out = vec.copy().astype(float)
    v = vec[keep]; cv = covar[keep]
    q = pd.qcut(pd.Series(cv).rank(method="first"), n_q, labels=False)
    med = pd.Series(v).groupby(q.values).transform("median").values
    glob = np.median(v[v > 0]) if (v > 0).any() else 1.0
    corr = np.where(med > 0, v / med * glob, 0.0)
    out[keep] = corr
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage2_dir", required=True)
    ap.add_argument("--samples", nargs="+", default=["520_520"])
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--outroot", required=True)
    ap.add_argument("--spatial_win", type=int, default=5)
    ap.add_argument("--normal_spots", default=None,
                    help="TSV (cols x_id,y_id) of reference/normal spots; the tumor aggregate is "
                         "divided by this normal aggregate (like the paper's /bulk-WGS step). "
                         "Use {sample} placeholder for per-sample files. If omitted, normalizes to genome median.")
    args = ap.parse_args()

    bins, off, n_bins = build_bins()
    gcf = pd.read_csv(os.path.join(args.ref_dir, "bin_gc_1mb.tsv"), sep="\t")
    gkey = (bins.chr + ":" + (bins.start // BIN).astype(str))
    gc = gkey.map(dict(zip(gcf.chr + ":" + (gcf.start // BIN).astype(str), gcf.gc))).fillna(0).values
    mapp = bed_frac_per_bin(os.path.join(args.ref_dir, "k100.umap.bed.gz"), off, n_bins, True) / BIN
    blk = bed_frac_per_bin(os.path.join(args.ref_dir, "hg38-blacklist.v2.bed.gz"), off, n_bins, False) / BIN

    for s in args.samples:
        sd = None
        for cand in (s, s + "_xy_map", s + "_out"):
            d = os.path.join(args.stage2_dir, cand)
            if os.path.isfile(os.path.join(d, "top50x50_fragments_with_cb.tsv.gz")):
                sd = d; break
        sd = sd or os.path.join(args.stage2_dir, s)
        frag = os.path.join(sd, "top50x50_fragments_with_cb.tsv.gz")
        if not os.path.isfile(frag):
            print("[%s] SKIP missing" % s); continue
        outdir = os.path.join(args.outroot, s); os.makedirs(outdir, exist_ok=True)
        counts = load_counts(frag, os.path.join(sd, "matched_spot_barcodes.tsv"), off, n_bins)
        tot = counts.sum(1).astype(float); present = tot > 0
        bulk = counts.sum(0).astype(float)

        npath = args.normal_spots.format(sample=s) if (args.normal_spots and "{sample}" in args.normal_spots) else args.normal_spots
        normal_mask = load_normal_mask(npath)

        # bin filter: mappable, not blacklisted, GC sane, adequate coverage
        cov_ok = bulk > np.median(bulk[bulk > 0]) * 0.25
        keep = (mapp > 0.5) & (blk < 0.3) & (gc > 0.3) & (gc < 0.65) & cov_ok
        print("[%s] bins kept=%d/%d" % (s, keep.sum(), n_bins))

        # ---- bulk CNA: GC then mappability quantile-correction ----
        b1 = quantile_correct(bulk, gc, keep)
        b2 = quantile_correct(b1, mapp, keep)
        relb = b2[keep] / np.median(b2[keep])
        # normal anchoring: divide the (corrected) tumor aggregate by the (corrected) normal-spot
        # aggregate -> removes per-bin technical pattern shared by both (paper's /bulk-WGS analogue).
        if normal_mask is not None and int((normal_mask & present).sum()) >= 20:
            sd_before = np.nanstd(np.log2(np.clip(relb, 1e-3, None)))
            nbulk = counts[normal_mask & present].sum(0).astype(float)
            nb2 = quantile_correct(quantile_correct(nbulk, gc, keep), mapp, keep)
            reln = nb2[keep] / np.median(nb2[keep]); reln[reln <= 0] = np.nan
            relb = relb / reln
            sd_after = np.nanstd(np.log2(np.clip(relb, 1e-3, None)))
            print("[%s] normal-anchored bulk CNA: %d reference spots (log2 SD %.3f -> %.3f)"
                  % (s, int((normal_mask & present).sum()), sd_before, sd_after))
            if sd_after < 0.4 * sd_before:
                print("[%s] *** SIGNAL-COLLAPSE WARNING: anchoring flattened the CNA (SD dropped to %.0f%%)."
                      % (s, 100 * sd_after / max(sd_before, 1e-9)))
                print("[%s] *** The reference shares the tumour's coverage pattern -> NO genuine internal normal."
                      % s)
                print("[%s] *** Keep the un-anchored profile; use an EXTERNAL reference (Stage-6 control / "
                      "histology-confirmed normal / RNA stroma) instead." % s)
        elif normal_mask is not None:
            print("[%s] WARN normal_spots too few; falling back to genome-median normalization" % s)
        rel = np.full(n_bins, np.nan)
        rel[keep] = relb
        log2 = np.log2(np.clip(rel, 1e-3, None))
        log2 = log2 - np.nanmedian(log2[keep])           # modal = diploid baseline
        gw = genome_smooth_1d(log2, bins, keep, 3)
        corr_gc = np.corrcoef(gw[keep], gc[keep])[0, 1]

        # ---- per-spot profiles (for spatial burden) ----
        C = counts[:, keep].astype(float)
        C = spatial_smooth(C, present, args.spatial_win)
        # per-spot GC + mappability quantile correction
        Cc = np.array([quantile_correct(C[i], gc[keep], np.ones(keep.sum(), bool)) if present[i] else C[i] for i in range(len(C))])
        Cc = np.array([quantile_correct(Cc[i], mapp[keep], np.ones(keep.sum(), bool)) if present[i] else Cc[i] for i in range(len(Cc))])
        t2 = Cc.sum(1); med = np.median(t2[present]); libn = np.zeros_like(Cc); libn[present] = Cc[present] / t2[present, None] * med
        ref = np.median(libn[present], 0); ref[ref <= 0] = np.nan
        spro = np.full_like(libn, np.nan); spro[present] = np.log2(np.clip(libn[present] / ref, 1e-3, None))
        spro = np.nan_to_num(spro) - np.nanmedian(np.nan_to_num(spro), axis=1, keepdims=True)

        # burden = correlation of each spot profile with the bulk tumor CNA profile
        bulk_vec = gw[keep]
        burden = np.full(GRID * GRID, np.nan)
        for i in np.where(present)[0]:
            v = spro[i]
            if np.std(v) > 0:
                burden[i] = np.corrcoef(v, bulk_vec)[0, 1]

        _save_outputs(s, outdir, bins, keep, gw, gc, mapp, burden, present, tot, corr_gc)
        print("[%s] DONE -> %s (bulk corr(CNA,GC)=%.2f after correction)" % (s, outdir, corr_gc))


def genome_smooth_1d(v, bins, keep, win):
    out = v.copy()
    for c in CHROMS:
        idx = np.where((bins.chr.values == c) & keep)[0]
        if len(idx) >= win:
            out[idx] = median_filter(v[idx], size=win, mode="nearest")
    return out


def spatial_smooth(mat, present, win):
    nb = mat.shape[1]
    data = np.zeros((GRID, GRID, nb)); mask = np.zeros((GRID, GRID))
    fl = np.nan_to_num(mat)
    for s in range(GRID * GRID):
        if present[s]:
            data[s // GRID, s % GRID] = fl[s]; mask[s // GRID, s % GRID] = 1
    num = uniform_filter(data * mask[:, :, None], size=(win, win, 1), mode="nearest")
    den = uniform_filter(mask, size=(win, win), mode="nearest")[:, :, None]
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0).reshape(GRID * GRID, nb)


def _save_outputs(s, outdir, bins, keep, gw, gc, mapp, burden, present, tot, corr_gc):
    kb = bins[keep].reset_index(drop=True).copy()
    kb["log2"] = gw[keep]
    kb.to_csv(os.path.join(outdir, "bulk_cna_profile.tsv"), sep="\t", index=False)
    # chrom boundaries
    bnds, centers, order, pos = [], [], [], 0
    for c in CHROMS:
        n = int((kb.chr.values == c).sum())
        if n:
            bnds.append(pos); centers.append(pos + n / 2); order.append(c.replace("chr", "")); pos += n
    bnds.append(pos)
    gene_x = {}
    for g, (gc_, mb) in PDAC_GENES.items():
        sel = kb.index[(kb.chr == gc_) & (kb.start <= mb * 1e6) & (kb.end > mb * 1e6)]
        if len(sel):
            gene_x[g] = sel[0]
    # bulk profile plot
    plt.figure(figsize=(15, 3.2))
    y = kb["log2"].values
    plt.scatter(range(len(y)), y, s=3, c=np.where(y > 0, "#c0392b", "#2471a3"))
    plt.axhline(0, color="grey", lw=0.6); plt.ylim(-1, 1)
    for b in bnds:
        plt.axvline(b - 0.5, color="k", lw=0.3, alpha=0.4)
    for g, gx in gene_x.items():
        plt.axvline(gx, color="darkorange", lw=0.7, ls="--", alpha=0.8)
        plt.text(gx, 0.92, g, fontsize=8, rotation=90, va="top", color="darkorange")
    plt.xticks(centers, order, fontsize=8); plt.ylabel("log2 (modal=0)")
    plt.title("%s: BULK tumor CNA (GC+mappability corrected, blacklist removed)" % s)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "bulk_cna_profile.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "bulk_cna_profile.pdf")); plt.close()
    # spatial burden map
    grid = np.full((GRID, GRID), np.nan)
    for sidx in range(GRID * GRID):
        grid[sidx // GRID, sidx % GRID] = burden[sidx]
    plt.figure(figsize=(7.5, 6))
    im = plt.imshow(grid, origin="lower", cmap="RdYlBu_r", interpolation="nearest", vmin=-0.5, vmax=0.8)
    plt.colorbar(im, label="similarity to bulk-tumor CNA  (tumor-content proxy)")
    plt.xlabel("Y"); plt.ylabel("X")
    plt.title("%s: spatial CNA-burden / tumor-likeness map" % s)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "spatial_cna_burden_map.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "spatial_cna_burden_map.pdf")); plt.close()
    pd.DataFrame({"x_id": np.repeat(np.arange(1, GRID + 1), GRID),
                 "y_id": np.tile(np.arange(1, GRID + 1), GRID),
                 "total_frags": tot.astype(int), "cna_burden": burden}).to_csv(
        os.path.join(outdir, "spot_cna_burden.tsv"), sep="\t", index=False)
    # per-chrom summary
    with open(os.path.join(outdir, "bulk_chrom_summary.tsv"), "w") as f:
        f.write("chr\tmean_log2\n")
        for c in CHROMS:
            v = kb[kb.chr == c]["log2"].mean()
            f.write("%s\t%.4f\n" % (c, v if not np.isnan(v) else 0))
        f.write("#corr_log2_vs_GC\t%.3f\n" % corr_gc)


if __name__ == "__main__":
    main()
