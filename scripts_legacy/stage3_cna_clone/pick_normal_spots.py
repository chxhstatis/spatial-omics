#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pick_normal_spots.py

Data-driven pseudo-normal selection (no histology needed): pick the FLATTEST spots
— lowest CNV burden, i.e. closest to diploid everywhere — as the reference set that
anchors the pseudo-normal in stage3_bulk_cna.py / stage3_per_spot_cnv.py.

Rationale: in a tumour with normal/stromal regions, those spots are near-diploid
(low |CN-2| fraction), so the flattest spots approximate a normal reference.

!!! VALIDITY CAVEAT (important) !!!
"Flattest" is defined RELATIVE TO the all-spot median. If the tumour is spread
uniformly across the section at low purity, that median is itself tumour-contaminated,
so the flattest spots = the ones most like the (tumour-bearing) average -> anchoring
to them CANCELS the real CNA (signal collapses to ~0). This auto-pick is only valid
when a genuine near-diploid subpopulation exists (e.g. a section containing adjacent
normal tissue). Pass --cnv_matrix to enable a self-check that DETECTS the circular
case and warns. For a uniformly low-purity section the correct reference is EXTERNAL:
a Stage-6 control sample (Panel-of-Normals), a histology-confirmed normal region (only
if one exists), or RNA-defined stroma (Stage-4). The --normal_spots interface in
stage3_bulk_cna.py / stage3_per_spot_cnv.py is built to take any of those.

Input : per_spot_cnv_burden.tsv  (cols x_id, y_id, total_frags, cnv_burden)
        [optional] --cnv_matrix per_spot_cnv_matrix.tsv.gz  -> enables the self-check
Output: normal_spots.tsv         (cols x_id, y_id) — the selected reference spots
        normal_spots_map.png     — where the reference spots sit on the 50x50 grid
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GRID = 50


def main():
    ap = argparse.ArgumentParser(description="Pick flattest spots as data-driven pseudo-normal reference.")
    ap.add_argument("--burden", required=True, help="per_spot_cnv_burden.tsv from stage3_per_spot_cnv.py")
    ap.add_argument("--out", required=True, help="output dir (writes normal_spots.tsv + map)")
    ap.add_argument("--quantile", type=float, default=0.30,
                    help="fraction of eligible spots (lowest CNV burden) to use as reference (default 0.30)")
    ap.add_argument("--min_frags", type=int, default=200,
                    help="only consider spots with at least this many fragments (stable reference)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    df = pd.read_csv(a.burden, sep="\t").dropna(subset=["cnv_burden"])
    elig = df[df.total_frags >= a.min_frags].copy()
    if len(elig) < 20:
        # relax the coverage floor if too few spots qualify
        elig = df.sort_values("total_frags", ascending=False).head(max(50, int(0.2 * len(df)))).copy()
        print("[pick_normal] few high-coverage spots; relaxed to top-coverage %d spots" % len(elig))
    thr = elig.cnv_burden.quantile(a.quantile)
    sel = elig[elig.cnv_burden <= thr][["x_id", "y_id"]].copy()
    sel.to_csv(os.path.join(a.out, "normal_spots.tsv"), sep="\t", index=False)
    print("[pick_normal] eligible=%d  selected=%d reference spots (cnv_burden<=%.3f, q=%.2f, min_frags=%d)"
          % (len(elig), len(sel), thr, a.quantile, a.min_frags))
    print("[pick_normal] NOTE: validity depends on a genuine near-diploid subpopulation existing. "
          "stage3_bulk_cna.py reports a 'signal-collapse' warning if anchoring to this set flattens the "
          "real CNA (the symptom of no internal normal -> use an external control). See docstring.")

    grid = np.full((GRID, GRID), np.nan)
    for r in df.itertuples():
        grid[int(r.x_id) - 1, int(r.y_id) - 1] = r.cnv_burden
    selmask = np.zeros((GRID, GRID))
    for r in sel.itertuples():
        selmask[int(r.x_id) - 1, int(r.y_id) - 1] = 1
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    im0 = ax[0].imshow(grid, origin="lower", cmap="viridis", interpolation="nearest")
    plt.colorbar(im0, ax=ax[0], label="CNV burden"); ax[0].set_title("per-spot CNV burden"); ax[0].set_xlabel("Y"); ax[0].set_ylabel("X")
    ax[1].imshow(selmask, origin="lower", cmap="Greys", interpolation="nearest")
    ax[1].set_title("selected pseudo-normal reference spots (n=%d)" % len(sel)); ax[1].set_xlabel("Y"); ax[1].set_ylabel("X")
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "normal_spots_map.png"), dpi=200); plt.close()
    print("[pick_normal] DONE -> %s/normal_spots.tsv" % a.out)


if __name__ == "__main__":
    main()
