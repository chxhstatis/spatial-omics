#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_bulk_atac_qc.py  (simple reconstruction)

Run inside a  <sample>_out/  directory. Reads the QC files produced by
stage1_preprocess_qc.sh and writes a QC summary figure + a text summary:
  - fragment length distribution  (tlen.hist.q30.txt)
  - TSS enrichment profile         (tss_profile_10bp.q30.tsv)   [expected ~flat for genomic DNA]
  - scalar metrics                 (flagstat*, frip.q30.txt, tss_enrichment.q30.txt, barcode_drop_stats.tsv)

Missing inputs are skipped gracefully.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_two_col(path, sep=None):
    xs, ys = [], []
    if not os.path.isfile(path):
        return xs, ys
    with open(path) as f:
        for ln in f:
            p = ln.split() if sep is None else ln.rstrip("\n").split(sep)
            if len(p) < 2:
                continue
            try:
                a, b = float(p[0]), float(p[1])
            except ValueError:
                continue
            xs.append(a)
            ys.append(b)
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--out", default="qc_summary.png")
    args = ap.parse_args()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # fragment length: file is "count tlen" (uniq -c output)
    cnt, tlen = read_two_col("tlen.hist.q30.txt")  # col1=count, col2=tlen
    if tlen:
        axes[0].plot(tlen, cnt, lw=0.8)
        if args.logy:
            axes[0].set_yscale("log")
        axes[0].set_xlabel("Fragment length (bp)")
        axes[0].set_ylabel("Count")
        axes[0].set_title("Fragment length distribution")
    else:
        axes[0].set_title("Fragment length (no data)")

    # TSS profile: header "offset\tcount"
    off, c = read_two_col("tss_profile_10bp.q30.tsv", sep="\t")
    if off:
        axes[1].plot(off, c, lw=0.8)
        axes[1].axvline(0, ls="--", alpha=0.3)
        axes[1].set_xlabel("Distance to TSS (bp)")
        axes[1].set_ylabel("Cut sites")
        axes[1].set_title("TSS profile (genomic DNA: expected ~flat)")
    else:
        axes[1].set_title("TSS profile (no data)")

    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    plt.savefig(os.path.splitext(args.out)[0] + ".pdf")
    plt.close()

    # text summary
    lines = []
    for path, label in [("tss_enrichment.q30.txt", None),
                        ("frip.q30.txt", None),
                        ("barcode_drop_stats.tsv", "barcode")]:
        if os.path.isfile(path):
            with open(path) as f:
                lines.append("# " + path)
                lines.extend(x.rstrip("\n") for x in f)
                lines.append("")
    with open("qc_summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("[plot_bulk_atac_qc] wrote", args.out, "and qc_summary.txt")


if __name__ == "__main__":
    main()
