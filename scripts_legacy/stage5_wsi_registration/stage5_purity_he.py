#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage5_purity_he.py

Independent, IMAGE-DERIVED tumour-cellularity / stromal-content estimate from the
matched H&E — orthogonal to the sequencing data. This BREAKS THE CIRCULARITY of
inferring "low purity" from low CNA amplitude and then using low purity to explain it:
here purity context is measured from morphology, not from the copy-number signal.

Method (within the drawn ROI box only):
  1. detect the ROI rectangle (red/blue line) -> mask the captured region.
  2. colour-deconvolve H&E (skimage rgb2hed): Haematoxylin (nuclei) + Eosin (stroma/cytoplasm).
  3. tissue mask = not white slide, not black scanner border.
  4. nuclei mask = haematoxylin above an Otsu threshold (within tissue).
  5. report:
       nuclei_fraction  = nuclei area / tissue area            (cellularity proxy; higher ~ more cellular/tumour)
       stroma_fraction  = eosin-high & haematoxylin-low area   (desmoplastic stroma proxy; higher ~ lower purity)
       mean_H, mean_E, H_over_E

CAVEAT (state in any report): this is a coarse cellularity/stromal-content proxy, NOT a
validated tumour-purity caller — stromal fibroblasts also contribute nuclei. In PDAC,
however, desmoplastic stroma is hypocellular dense collagen (eosin-rich, nuclei-poor),
so low nuclei_fraction / high stroma_fraction independently indicates a stroma-dominated,
low-tumour-content section. For a calibrated purity, pair with pathologist cellularity or
RNA deconvolution (RCTD). The value here is providing an INDEPENDENT axis to test the
purity hypothesis against, not a final purity number.

Outputs (per sample): purity_he.tsv, he_deconv_panel.png
Cohort: purity_he_cohort.tsv
"""
import argparse, os
import numpy as np, pandas as pd
import cv2
from skimage.color import rgb2hed
from skimage.filters import threshold_otsu
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def detect_box_mask(img):
    """largest red/blue/green drawn rectangle -> filled polygon mask + colour."""
    b, g, r = img[:, :, 0].astype(int), img[:, :, 1].astype(int), img[:, :, 2].astype(int)
    masks = {"red": (r > 120) & (g < 90) & (b < 90),
             "blue": (b > 110) & (r < 90) & (g < 90),
             "green": (g > 110) & (r < 90) & (b < 90)}
    color = max(masks, key=lambda k: int(masks[k].sum()))
    m = (masks[color].astype(np.uint8)) * 255
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, np.ones((15, 15), np.uint8))
    cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    c = max(cnts, key=cv2.contourArea)
    box = cv2.boxPoints(cv2.minAreaRect(c)).astype(np.int32)
    poly = np.zeros(img.shape[:2], np.uint8)
    cv2.fillPoly(poly, [box], 1)
    return poly.astype(bool), color, box


def main():
    ap = argparse.ArgumentParser(description="H&E-derived cellularity/stroma proxy (independent purity anchor).")
    ap.add_argument("--he", required=True)
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    img = cv2.imread(a.he)
    if img is None:
        raise RuntimeError("cannot read %s" % a.he)
    roi, color, box = detect_box_mask(img)
    rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB).astype(float) / 255.0

    # colour deconvolution
    hed = rgb2hed(rgb)
    H = hed[:, :, 0]; E = hed[:, :, 1]
    # normalise H,E to [0,1] within ROI for thresholding/visual
    def norm(x, mask):
        lo, hi = np.percentile(x[mask], 1), np.percentile(x[mask], 99)
        return np.clip((x - lo) / (hi - lo + 1e-9), 0, 1)
    Hn = norm(H, roi); En = norm(E, roi)

    mx = rgb.max(2); mn = rgb.min(2)
    white = mn > 0.80
    black = mx < 0.20
    tissue = roi & ~white & ~black
    n_tissue = int(tissue.sum())
    if n_tissue < 1000:
        raise RuntimeError("too little tissue inside ROI")

    # nuclei = haematoxylin above Otsu (within tissue)
    thr = threshold_otsu(Hn[tissue])
    nuclei = tissue & (Hn > thr)
    stroma = tissue & (En > np.percentile(En[tissue], 60)) & (Hn < thr)

    nuclei_fraction = float(nuclei.sum()) / n_tissue
    stroma_fraction = float(stroma.sum()) / n_tissue
    mean_H = float(np.mean(Hn[tissue])); mean_E = float(np.mean(En[tissue]))
    tissue_fraction_of_roi = n_tissue / float(roi.sum())

    row = dict(sample=a.sample, box_color=color,
               tissue_frac_of_ROI=round(tissue_fraction_of_roi, 3),
               nuclei_fraction=round(nuclei_fraction, 3),
               stroma_fraction=round(stroma_fraction, 3),
               mean_haematoxylin=round(mean_H, 3), mean_eosin=round(mean_E, 3),
               H_over_E=round(mean_H / (mean_E + 1e-9), 3))
    pd.DataFrame([row]).to_csv(os.path.join(a.out, "purity_he.tsv"), sep="\t", index=False)

    # panel figure
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))
    show = rgb.copy()
    cv2.polylines((show * 255).astype(np.uint8), [box], True, (255, 255, 0), 3)
    ax[0].imshow(rgb); ax[0].plot(*np.vstack([box, box[0]]).T, color="yellow", lw=1.5); ax[0].set_title("%s H&E + ROI" % a.sample)
    ax[1].imshow(np.where(tissue, Hn, np.nan), cmap="Purples"); ax[1].set_title("Haematoxylin (nuclei) in ROI")
    ax[2].imshow(np.where(tissue, En, np.nan), cmap="Reds"); ax[2].set_title("Eosin (stroma) in ROI")
    overlay = np.zeros((*img.shape[:2], 3))
    overlay[nuclei] = [0.4, 0.0, 0.6]; overlay[stroma] = [0.95, 0.6, 0.6]
    ax[3].imshow(overlay); ax[3].set_title("nuclei(purple)=%.0f%%  stroma(pink)=%.0f%%" % (100*nuclei_fraction, 100*stroma_fraction))
    for x in ax: x.axis("off")
    plt.suptitle("%s — H&E-derived cellularity (independent of sequencing). nuclei_frac=%.2f  stroma_frac=%.2f  H/E=%.2f"
                 % (a.sample, nuclei_fraction, stroma_fraction, row["H_over_E"]))
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "he_deconv_panel.png"), dpi=130); plt.close()
    print("[%s] nuclei_frac=%.3f stroma_frac=%.3f mean_H=%.3f mean_E=%.3f H/E=%.3f -> %s"
          % (a.sample, nuclei_fraction, stroma_fraction, mean_H, mean_E, row["H_over_E"], a.out))


if __name__ == "__main__":
    main()
