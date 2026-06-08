#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_samples_cna.py

Cohort-level comparison of bulk CNA across multiple tumour samples (no normal needed).
Reads each sample's bulk_cna_profile.tsv (chr,start,end,log2) from a stage3_out_bulk
output root, aligns bins, and produces:

  cohort_cna_profiles.png   stacked per-sample CNA profiles aligned along the genome
  recurrence_track.png      per-bin recurrence: how many samples gain / lose each bin
  sample_correlation.png    sample x sample CNA-profile correlation (which tumours are alike)
  recurrence.tsv            chr,start,end, n_gain, n_loss, mean_log2  (recurrent events table)
  sample_correlation.tsv    pairwise correlation matrix
  recurrent_genes.tsv       canonical PDAC loci: per-sample log2 + #gain/#loss across cohort

Use to find which CNA events recur across patients (PDAC drivers) vs private events.
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
LOCI = {"CDKN2A(9p)": ("chr9", 21_000_000), "MYC(8q)": ("chr8", 127_000_000),
        "KRAS(12p)": ("chr12", 25_000_000), "TP53(17p)": ("chr17", 7_000_000),
        "SMAD4(18q)": ("chr18", 51_000_000), "GATA6(18q11)": ("chr18", 22_000_000)}


def main():
    ap = argparse.ArgumentParser(description="Cross-sample (cohort) bulk CNA comparison.")
    ap.add_argument("--bulk_root", required=True, help="stage3_out_bulk dir (holds <sample>/bulk_cna_profile.tsv)")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--thresh", type=float, default=0.2, help="|log2| gain/loss call threshold")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    # load + align on (chr,start)
    mats = {}
    for s in a.samples:
        p = os.path.join(a.bulk_root, s, "bulk_cna_profile.tsv")
        if not os.path.isfile(p):
            print("[compare] SKIP %s (no %s)" % (s, p)); continue
        d = pd.read_csv(p, sep="\t")
        mats[s] = d.set_index(d.chr + ":" + d.start.astype(str))["log2"]
    samples = [s for s in a.samples if s in mats]
    ref = pd.read_csv(os.path.join(a.bulk_root, samples[0], "bulk_cna_profile.tsv"), sep="\t")
    key = ref.chr + ":" + ref.start.astype(str)
    M = pd.DataFrame({s: mats[s].reindex(key).values for s in samples})
    M.insert(0, "chr", ref.chr.values); M.insert(1, "start", ref.start.values); M.insert(2, "end", ref.end.values)
    # genome order
    order = {c: i for i, c in enumerate(CHROMS)}
    M = M[M.chr.isin(CHROMS)].copy()
    M["_o"] = M.chr.map(order); M = M.sort_values(["_o", "start"]).reset_index(drop=True)
    L = M[samples].values  # bins x samples

    # chrom boundaries for plotting
    bnds, centers, labs, pos = [], [], [], 0
    for c in CHROMS:
        m = int((M.chr == c).sum())
        if m: bnds.append(pos); centers.append(pos + m / 2); labs.append(c.replace("chr", "")); pos += m
    def gene_x(c, p):
        sub = M[(M.chr == c) & (M.start <= p) & (M.end > p)]
        return sub.index[0] if len(sub) else None

    # ---- fig 1: stacked profiles ----
    fig, axes = plt.subplots(len(samples), 1, figsize=(15, 1.8 * len(samples)), squeeze=False, sharex=True)
    for ax, s in zip(axes[:, 0], samples):
        v = M[s].values
        ax.scatter(np.arange(len(v)), v, s=2, c=np.where(v > a.thresh, "#c0392b", np.where(v < -a.thresh, "#2471a3", "#b0b0b0")))
        ax.axhline(0, color="grey", lw=0.5); ax.set_ylim(-1, 1)
        for b in bnds: ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.3)
        ax.set_ylabel(s, fontsize=9)
    for k, (c, p) in LOCI.items():
        x = gene_x(c, p)
        if x is not None:
            for ax in axes[:, 0]: ax.axvline(x, color="orange", lw=0.6, ls=":", alpha=0.6)
            axes[0, 0].text(x, 1.02, k, rotation=90, fontsize=6, color="darkorange", va="bottom")
    axes[-1, 0].set_xticks(centers); axes[-1, 0].set_xticklabels(labs, fontsize=7)
    axes[0, 0].set_title("Cohort bulk CNA — per-sample log2 (red=gain, blue=loss)")
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "cohort_cna_profiles.png"), dpi=200); plt.close()

    # ---- fig 2: recurrence track ----
    n_gain = np.nansum(L > a.thresh, axis=1); n_loss = np.nansum(L < -a.thresh, axis=1)
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.fill_between(np.arange(len(n_gain)), 0, n_gain, color="#c0392b", step="mid", label="gain")
    ax.fill_between(np.arange(len(n_loss)), 0, -n_loss, color="#2471a3", step="mid", label="loss")
    for b in bnds: ax.axvline(b - 0.5, color="k", lw=0.3, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5); ax.set_ylim(-len(samples), len(samples))
    ax.set_yticks(range(-len(samples), len(samples) + 1))
    ax.set_yticklabels([str(abs(t)) for t in range(-len(samples), len(samples) + 1)])
    ax.set_ylabel("# samples (loss / gain)"); ax.set_xticks(centers); ax.set_xticklabels(labs, fontsize=7)
    for k, (c, p) in LOCI.items():
        x = gene_x(c, p)
        if x is not None: ax.axvline(x, color="orange", lw=0.6, ls=":", alpha=0.7); ax.text(x, len(samples)*0.95, k, rotation=90, fontsize=6, color="darkorange", va="top")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("CNA recurrence across %d tumours (|log2|>%.2f)" % (len(samples), a.thresh))
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "recurrence_track.png"), dpi=200); plt.close()

    # ---- fig 3: sample correlation (AUTOSOMES ONLY: chrX would cluster samples by SEX, not tumour biology) ----
    auto_mask = (M.chr != "chrX").values
    C = np.corrcoef(np.nan_to_num(L[auto_mask]).T)
    chrx_mean = {s: float(np.nanmean(M.loc[M.chr == "chrX", s].values)) for s in samples}
    sex = {s: ("male" if chrx_mean[s] < -0.3 else "female") for s in samples}
    fig, ax = plt.subplots(figsize=(1.2 * len(samples) + 2, 1.2 * len(samples) + 1))
    im = ax.imshow(C, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1))
    ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(samples))); ax.set_yticklabels(samples, fontsize=8)
    for i in range(len(samples)):
        for j in range(len(samples)): ax.text(j, i, "%.2f" % C[i, j], ha="center", va="center", fontsize=8)
    plt.colorbar(im, label="CNA-profile correlation"); ax.set_title("Sample similarity (autosomal bulk CNA; chrX excluded)")
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "sample_correlation.png"), dpi=200); plt.close()
    pd.DataFrame({"sample": samples, "chrX_mean_log2": [round(chrx_mean[s], 2) for s in samples],
                  "inferred_sex": [sex[s] for s in samples]}).to_csv(
        os.path.join(a.out, "inferred_sex.tsv"), sep="\t", index=False)

    # ---- tables ----
    rec = M[["chr", "start", "end"]].copy()
    rec["n_gain"] = n_gain; rec["n_loss"] = n_loss; rec["mean_log2"] = np.nanmean(L, axis=1)
    rec.to_csv(os.path.join(a.out, "recurrence.tsv"), sep="\t", index=False)
    pd.DataFrame(C, index=samples, columns=samples).to_csv(os.path.join(a.out, "sample_correlation.tsv"), sep="\t")
    rows = []
    for k, (c, p) in LOCI.items():
        x = gene_x(c, p)
        if x is None: continue
        vals = L[x]
        rows.append(dict(locus=k, **{s: round(float(M[s].values[x]), 3) for s in samples},
                         n_gain=int(np.nansum(vals > a.thresh)), n_loss=int(np.nansum(vals < -a.thresh))))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "recurrent_genes.tsv"), sep="\t", index=False)

    print("[compare] %d samples, %d aligned bins -> %s" % (len(samples), len(M), a.out))
    print("[compare] inferred sex (chrX):", {s: sex[s] for s in samples})
    print("[compare] AUTOSOMAL sample correlation (chrX excluded):\n",
          pd.DataFrame(C, index=samples, columns=samples).round(2).to_string())


if __name__ == "__main__":
    main()
