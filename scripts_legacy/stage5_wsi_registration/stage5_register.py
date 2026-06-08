#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage5_register.py

Register an H&E overview image (chip-in-place, with a drawn rectangular ROI marking
the 50x50 capture region) to the DBiT spot grid, and read per-spot tissue density.

Pipeline:
  1. Detect the drawn ROI rectangle (red or blue line) -> 4 corners (cv2.minAreaRect,
     handles rotated boxes such as sample A).
  2. The grid->image mapping has an 8-fold orientation ambiguity (4 rotations x 2 flips)
     that no fiducial resolves. We resolve it FROM THE DATA: the correct orientation is
     the one whose per-spot H&E tissue density best correlates with per-spot DNA coverage
     (more tissue -> more fragments). We try all 8 and pick the max-correlation mapping.
     The chosen correlation r is the registration's validation metric.
  3. Emit: the perspective transform, per-spot pixel coords + tissue density, an H&E+grid
     overlay for visual QC, and a density-vs-coverage validation figure.

Inputs:
  --he <image.jpg>             H&E overview (the .jpg next to the .kfb)
  --coverage <tsv>             per-spot coverage with cols x_id,y_id,total_frags
                               (e.g. stage3_out_perspot/<sample>/per_spot_cnv_burden.tsv)
Outputs (in --out):
  spot_to_image_perspective.npy   3x3 matrix mapping grid (u,v in [0,1]) -> image px
  spot_pixel_coords.tsv           x_id,y_id,px,py,tissue_density
  spot_density.tsv                x_id,y_id,tissue_density  (for downstream --normal_spots etc.)
  chip_grid_overlay.png           H&E with the 50x50 grid + ROI drawn
  density_vs_coverage.png         side-by-side grid maps + scatter (validation)
  register_summary.tsv            box color, corners, chosen orientation, validation r
"""
import argparse, os
import numpy as np, pandas as pd
import cv2
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GRID = 50


def detect_box(img):
    """return (corners 4x2 float, color) of the drawn ROI rectangle (red or blue line)."""
    b, g, r = img[:, :, 0].astype(int), img[:, :, 1].astype(int), img[:, :, 2].astype(int)
    masks = {"red": (r > 120) & (g < 90) & (b < 90),
             "blue": (b > 110) & (r < 90) & (g < 90),
             "green": (g > 110) & (r < 90) & (b < 90)}
    color = max(masks, key=lambda k: int(masks[k].sum()))
    m = masks[color].astype(np.uint8) * 255
    if m.sum() < 50 * 255:
        raise RuntimeError("no colored ROI box found (red/blue/green)")
    # close gaps along the drawn line, then take the LARGEST contour (the ROI outline) only,
    # so scattered specks / labels of the same colour don't distort the rotated rect.
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, np.ones((15, 15), np.uint8))
    cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        raise RuntimeError("no ROI contour")
    c = max(cnts, key=cv2.contourArea)
    rect = cv2.minAreaRect(c)
    box = cv2.boxPoints(rect)  # 4x2, cyclic order
    return box.astype(np.float64), color


def orientations(box):
    """8 dihedral mappings of the unit square corners -> the 4 box corners."""
    UC = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], np.float64)  # (u,v)
    outs = []
    for rev in (box, box[::-1]):
        for k in range(4):
            outs.append((UC, np.roll(rev, k, axis=0)))
    return outs


def tissue_density(img, px, py, half):
    """fraction of a (2half x 2half) tile around (px,py) that is stained tissue
    (not white background, not black scanner border)."""
    H, W = img.shape[:2]
    x0, x1 = max(0, int(px - half)), min(W, int(px + half) + 1)
    y0, y1 = max(0, int(py - half)), min(H, int(py + half) + 1)
    if x1 <= x0 or y1 <= y0:
        return np.nan
    tile = img[y0:y1, x0:x1].astype(int)
    mx = tile.max(2); mn = tile.min(2)
    white = mn > 200            # background slide
    black = mx < 50             # scanner border / mounting
    tissue = ~white & ~black
    return float(tissue.mean())


def main():
    ap = argparse.ArgumentParser(description="Register H&E (ROI box) to the 50x50 DBiT grid.")
    ap.add_argument("--he", required=True)
    ap.add_argument("--coverage", required=True, help="tsv with x_id,y_id,total_frags")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    img = cv2.imread(a.he)
    if img is None:
        raise RuntimeError("cannot read %s" % a.he)
    box, color = detect_box(img)
    print("[%s] ROI box color=%s corners=%s" % (a.sample, color, np.round(box).tolist()))

    cov = pd.read_csv(a.coverage, sep="\t")
    cov = cov[cov.total_frags > 0].copy()
    u = (cov.x_id.values - 0.5) / GRID
    v = (cov.y_id.values - 0.5) / GRID
    cov_vec = cov.total_frags.values.astype(float)
    cov_log = np.log1p(cov_vec)
    pitch = np.linalg.norm(box[0] - box[1]) / GRID
    half = max(2.0, pitch / 2.0)

    # DBiT per-spot coverage is dominated by microfluidic channel-efficiency STRIPES, not tissue,
    # so we resolve orientation on the stripe-suppressed COARSE TISSUE ENVELOPE: smooth both the
    # H&E-density grid and the coverage grid (kills stripes + per-spot noise) and correlate those.
    from scipy.ndimage import median_filter
    def grid_of(vals, xi, yi):
        g = np.full((GRID, GRID), np.nan); g[xi - 1, yi - 1] = vals; return g
    def smooth(g):
        m = np.isfinite(g); gg = np.where(m, g, np.nanmean(g[m]))
        return median_filter(gg, size=5, mode="nearest"), m
    xi = cov.x_id.values.astype(int); yi = cov.y_id.values.astype(int)
    cov_sm, _ = smooth(grid_of(cov_log, xi, yi))

    best = None; all_r = []
    for oi, (UC, Q) in enumerate(orientations(box)):
        Hmat = cv2.getPerspectiveTransform(UC.astype(np.float32), Q.astype(np.float32))
        uv = np.column_stack([u, v, np.ones_like(u)])
        pix = (Hmat @ uv.T).T; pix = pix[:, :2] / pix[:, 2:3]
        dens = np.array([tissue_density(img, pix[i, 0], pix[i, 1], half) for i in range(len(pix))])
        if np.nanstd(dens) == 0:
            all_r.append(np.nan); continue
        dens_sm, m = smooth(grid_of(dens, xi, yi))
        r = np.corrcoef(dens_sm[m].ravel(), cov_sm[m].ravel())[0, 1]
        all_r.append(round(float(r), 3))
        if best is None or r > best["r"]:
            best = dict(oi=oi, H=Hmat, pix=pix, dens=dens, r=float(r))
    if best is None:
        raise RuntimeError("orientation search failed (no density variation)")
    rs = sorted([x for x in all_r if np.isfinite(x)], reverse=True)
    margin = (rs[0] - rs[1]) if len(rs) > 1 else np.nan
    resolved = bool(best["r"] > 0.25 and np.isfinite(margin) and margin > 0.1)
    print("[%s] orientation r over 8: %s" % (a.sample, all_r))
    print("[%s] chosen %d/8 r=%.3f margin=%.3f -> orientation %s"
          % (a.sample, best["oi"], best["r"], margin if np.isfinite(margin) else -1,
             "RESOLVED" if resolved else "AMBIGUOUS (geometry ok; spot-index orientation needs a fiducial/user)"))

    # ---- outputs ----
    np.save(os.path.join(a.out, "spot_to_image_perspective.npy"), best["H"])
    out = cov[["x_id", "y_id"]].copy()
    out["px"] = best["pix"][:, 0]; out["py"] = best["pix"][:, 1]; out["tissue_density"] = best["dens"]
    out.to_csv(os.path.join(a.out, "spot_pixel_coords.tsv"), sep="\t", index=False)
    out[["x_id", "y_id", "tissue_density"]].to_csv(os.path.join(a.out, "spot_density.tsv"), sep="\t", index=False)
    pd.DataFrame([dict(sample=a.sample, box_color=color, orientation=best["oi"],
                       validation_r=round(best["r"], 3), orientation_margin=round(margin, 3) if np.isfinite(margin) else np.nan,
                       orientation_resolved=resolved, all_orientation_r=str(all_r),
                       pitch_px=round(pitch, 1), corners=str(np.round(box).tolist()))]).to_csv(
        os.path.join(a.out, "register_summary.tsv"), sep="\t", index=False)

    # overlay: H&E (RGB) + grid lines + ROI box
    rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    plt.figure(figsize=(11, 11 * img.shape[0] / img.shape[1]))
    plt.imshow(rgb)
    Hm = best["H"]
    def to_px(uu, vv):
        p = Hm @ np.array([uu, vv, 1.0]); return p[:2] / p[2]
    for gi in range(0, GRID + 1, 5):
        a1 = to_px(gi / GRID, 0); a2 = to_px(gi / GRID, 1)
        b1 = to_px(0, gi / GRID); b2 = to_px(1, gi / GRID)
        plt.plot([a1[0], a2[0]], [a1[1], a2[1]], color="lime", lw=0.5, alpha=0.7)
        plt.plot([b1[0], b2[0]], [b1[1], b2[1]], color="lime", lw=0.5, alpha=0.7)
    bx = np.vstack([box, box[0]])
    plt.plot(bx[:, 0], bx[:, 1], color="yellow", lw=1.5)
    plt.title("%s: H&E + DBiT 50x50 grid (orientation %d, r=%.2f)" % (a.sample, best["oi"], best["r"]))
    plt.axis("off"); plt.tight_layout()
    plt.savefig(os.path.join(a.out, "chip_grid_overlay.png"), dpi=150); plt.close()

    # validation: density grid vs coverage grid + scatter
    dg = np.full((GRID, GRID), np.nan); cg = np.full((GRID, GRID), np.nan)
    for i in range(len(out)):
        dg[int(out.x_id.values[i]) - 1, int(out.y_id.values[i]) - 1] = best["dens"][i]
        cg[int(out.x_id.values[i]) - 1, int(out.y_id.values[i]) - 1] = cov_log[i]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    im0 = ax[0].imshow(dg, origin="lower", cmap="pink_r"); ax[0].set_title("H&E tissue density"); plt.colorbar(im0, ax=ax[0])
    im1 = ax[1].imshow(cg, origin="lower", cmap="viridis"); ax[1].set_title("DNA coverage (log frags)"); plt.colorbar(im1, ax=ax[1])
    ax[2].scatter(best["dens"], cov_log, s=4, alpha=0.4); ax[2].set_xlabel("H&E density"); ax[2].set_ylabel("log coverage")
    ax[2].set_title("validation r=%.3f" % best["r"])
    for x in (ax[0], ax[1]): x.set_xlabel("Y"); x.set_ylabel("X")
    plt.suptitle("%s: registration validation (density should track coverage)" % a.sample)
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "density_vs_coverage.png"), dpi=150); plt.close()
    print("[%s] DONE -> %s" % (a.sample, a.out))


if __name__ == "__main__":
    main()
