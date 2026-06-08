"""H&E ↔ spot-grid registration.

Map the spot grid onto a matched H&E overview that has a drawn rectangular ROI marking
the capture region. The geometry (ROI box → grid) is robust; the 8-fold spot-index
ORIENTATION (which corner is spot (1,1)) is resolved from data by correlating the
stripe-suppressed H&E tissue envelope with the DNA coverage envelope — but on DBiT data,
where coverage is dominated by microfluidic channel stripes and tissue fills the ROI, the
orientation is often AMBIGUOUS. The result flags this honestly (`orientation_resolved`);
a fiducial, the chip convention, or matched RNA resolves it.

Requires the optional ``he`` extra (opencv-python), imported lazily.
"""
from __future__ import annotations

import numpy as np


def _orientations(box):
    """8 dihedral mappings of the unit-square corners -> the 4 ROI-box corners."""
    UC = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], np.float64)
    return [(UC, np.roll(rev, k, axis=0)) for rev in (box, box[::-1]) for k in range(4)]


def _tissue_density(img, px, py, half):
    """fraction of the tile around (px,py) that is stained tissue (not white, not black)."""
    H, W = img.shape[:2]
    x0, x1 = max(0, int(px - half)), min(W, int(px + half) + 1)
    y0, y1 = max(0, int(py - half)), min(H, int(py + half) + 1)
    if x1 <= x0 or y1 <= y0:
        return np.nan
    tile = img[y0:y1, x0:x1].astype(int)
    mx, mn = tile.max(2), tile.min(2)
    return float((~(mn > 200) & ~(mx < 50)).mean())


def register(adata, image_path, *, coverage_key="total_frags", save_overlay=None):
    """Register a matched H&E (with a drawn ROI box) to the spot grid.

    Writes to ``adata``:
      - ``obsm['X_he']`` — per-spot (px, py) pixel coordinates in the H&E image;
      - ``obs['tissue_density']`` — H&E tissue fraction at each spot;
      - ``uns['spatial_omics_register']`` — transform (3×3), ROI box corners + colour,
        chosen orientation, ``validation_r``, ``orientation_resolved`` (+ all 8 scores).

    Parameters
    ----------
    image_path : str
        H&E overview with a red/blue/green rectangular ROI outline marking the 50×50 area.
    coverage_key : str
        ``adata.obs`` column used as DNA coverage for orientation scoring (default total_frags).
    save_overlay : str, optional
        If given, save an H&E + grid + ROI overlay PNG there (visual QC; orientation-robust).

    Returns the ``uns['spatial_omics_register']`` summary dict.
    """
    try:
        import cv2
    except ImportError as e:  # pragma: no cover
        raise ImportError("register needs the 'he' extra: pip install 'spatial-omics[he]'") from e
    from scipy.ndimage import median_filter

    from .he import _detect_box_mask

    img = cv2.imread(str(image_path))
    if img is None:
        raise FileNotFoundError(f"cannot read {image_path}")
    _, color, box = _detect_box_mask(img)
    box = box.astype(np.float64)

    xi = adata.obs["x_id"].values.astype(int)
    yi = adata.obs["y_id"].values.astype(int)
    nx, ny = xi.max(), yi.max()
    u = (xi - 0.5) / nx
    v = (yi - 0.5) / ny
    cov_log = np.log1p(adata.obs[coverage_key].values.astype(float))
    pitch = np.linalg.norm(box[0] - box[1]) / max(nx, ny)
    half = max(2.0, pitch / 2.0)

    def grid_of(vals):
        g = np.full((nx, ny), np.nan); g[xi - 1, yi - 1] = vals; return g

    def smooth(g):
        m = np.isfinite(g); gg = np.where(m, g, np.nanmean(g[m]))
        return median_filter(gg, size=5, mode="nearest"), m

    cov_sm, _ = smooth(grid_of(cov_log))

    best, all_r = None, []
    for oi, (UC, Q) in enumerate(_orientations(box)):
        Hmat = cv2.getPerspectiveTransform(UC.astype(np.float32), Q.astype(np.float32))
        uv = np.column_stack([u, v, np.ones_like(u)])
        pix = (Hmat @ uv.T).T
        pix = pix[:, :2] / pix[:, 2:3]
        dens = np.array([_tissue_density(img, pix[i, 0], pix[i, 1], half) for i in range(len(pix))])
        if np.nanstd(dens) == 0:
            all_r.append(np.nan); continue
        dens_sm, m = smooth(grid_of(dens))
        r = float(np.corrcoef(dens_sm[m].ravel(), cov_sm[m].ravel())[0, 1])
        all_r.append(round(r, 3))
        if best is None or r > best["r"]:
            best = {"oi": oi, "H": Hmat, "pix": pix, "dens": dens, "r": r}
    if best is None:
        raise RuntimeError("orientation search failed (no tissue-density variation in ROI)")

    rs = sorted([x for x in all_r if np.isfinite(x)], reverse=True)
    margin = (rs[0] - rs[1]) if len(rs) > 1 else np.nan
    resolved = bool(best["r"] > 0.25 and np.isfinite(margin) and margin > 0.1)

    adata.obsm["X_he"] = best["pix"]
    adata.obs["tissue_density"] = best["dens"]
    res = {"box_color": color, "box_corners": np.round(box).astype(int).tolist(),
           "transform": best["H"].tolist(), "orientation": int(best["oi"]),
           "validation_r": round(best["r"], 3),
           "orientation_margin": round(margin, 3) if np.isfinite(margin) else None,
           "orientation_resolved": resolved, "all_orientation_r": all_r,
           "pitch_px": round(float(pitch), 1),
           "note": ("orientation RESOLVED from data" if resolved else
                    "geometry OK but spot-index orientation AMBIGUOUS — needs a fiducial / chip "
                    "convention / matched RNA (DBiT coverage is channel-stripe-dominated)")}
    adata.uns["spatial_omics_register"] = res

    if save_overlay:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
        Hm = best["H"]
        topx = lambda uu, vv: (Hm @ np.array([uu, vv, 1.0]))[:2] / (Hm @ np.array([uu, vv, 1.0]))[2]
        plt.figure(figsize=(10, 10 * img.shape[0] / img.shape[1]))
        plt.imshow(rgb)
        for gi in range(0, max(nx, ny) + 1, 5):
            a1, a2 = topx(gi / nx, 0), topx(gi / nx, 1)
            b1, b2 = topx(0, gi / ny), topx(1, gi / ny)
            plt.plot([a1[0], a2[0]], [a1[1], a2[1]], color="lime", lw=0.5, alpha=0.7)
            plt.plot([b1[0], b2[0]], [b1[1], b2[1]], color="lime", lw=0.5, alpha=0.7)
        bx = np.vstack([box, box[0]])
        plt.plot(bx[:, 0], bx[:, 1], color="yellow", lw=1.5)
        plt.title("H&E + spot grid (orientation %d, r=%.2f%s)"
                  % (best["oi"], best["r"], "" if resolved else ", AMBIGUOUS"))
        plt.axis("off"); plt.tight_layout(); plt.savefig(save_overlay, dpi=150); plt.close()
    return res
