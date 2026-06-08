"""H&E morphology — an independent, image-derived tumour-cellularity / stromal-content
estimate. Orthogonal to the sequencing, it breaks the circularity of inferring "low
purity" from low CNA amplitude and then using purity to explain it.

Requires the optional ``he`` extra (``pip install 'spatial-omics[he]'`` ->
opencv-python + scikit-image); imported lazily so the core package has no such dependency.

Method (within the drawn ROI box only): colour-deconvolve H&E (skimage rgb2hed) into
Haematoxylin (nuclei) + Eosin (stroma/cytoplasm); nuclei = haematoxylin above an Otsu
threshold inside the tissue mask. Reports ``nuclei_fraction`` (cellularity proxy; higher
~ more cellular/tumour) and ``stroma_fraction`` (eosin-high & nuclei-low; desmoplastic
stroma proxy; higher ~ lower purity).

Caveat: a coarse cellularity/stromal-content proxy, NOT a calibrated tumour purity
(stromal fibroblasts also contribute nuclei; staining intensity varies between slides).
Use it as an independent axis to test the purity hypothesis, paired with the
``pp.normal_anchor`` signal-collapse guard and (when available) RNA deconvolution.
"""
from __future__ import annotations

import numpy as np


def _detect_box_mask(img):
    """Largest drawn red/blue/green rectangle -> filled polygon mask + colour + corners."""
    import cv2
    b, g, r = img[:, :, 0].astype(int), img[:, :, 1].astype(int), img[:, :, 2].astype(int)
    masks = {"red": (r > 120) & (g < 90) & (b < 90),
             "blue": (b > 110) & (r < 90) & (g < 90),
             "green": (g > 110) & (r < 90) & (b < 90)}
    color = max(masks, key=lambda k: int(masks[k].sum()))
    m = (masks[color].astype(np.uint8)) * 255
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, np.ones((15, 15), np.uint8))
    cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        raise RuntimeError("no colored ROI box found (red/blue/green)")
    c = max(cnts, key=cv2.contourArea)
    box = cv2.boxPoints(cv2.minAreaRect(c)).astype(np.int32)
    poly = np.zeros(img.shape[:2], np.uint8)
    cv2.fillPoly(poly, [box], 1)
    return poly.astype(bool), color, box


def he_purity(image_path, *, adata=None, eosin_pct=60):
    """Image-derived cellularity / stromal-content from a matched H&E (ROI box required).

    Parameters
    ----------
    image_path : str
        H&E overview image (jpg/png/tiff) with a drawn rectangular ROI marking the capture
        region (red/blue/green outline).
    adata : AnnData, optional
        If given, the result dict is also stored in ``adata.uns['spatial_omics_he_purity']``.
    eosin_pct : float
        Percentile of in-tissue eosin above which (with low haematoxylin) a pixel counts as
        stroma.

    Returns
    -------
    dict with ``nuclei_fraction``, ``stroma_fraction``, ``mean_haematoxylin``, ``mean_eosin``,
    ``H_over_E``, ``box_color``, ``tissue_frac_of_ROI``.
    """
    try:
        import cv2
        from skimage.color import rgb2hed
        from skimage.filters import threshold_otsu
    except ImportError as e:  # pragma: no cover
        raise ImportError("he_purity needs the 'he' extra: pip install 'spatial-omics[he]'") from e

    img = cv2.imread(str(image_path))
    if img is None:
        raise FileNotFoundError(f"cannot read {image_path}")
    roi, color, _ = _detect_box_mask(img)
    rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB).astype(float) / 255.0

    hed = rgb2hed(rgb)
    H, E = hed[:, :, 0], hed[:, :, 1]

    def _norm(x, mask):
        lo, hi = np.percentile(x[mask], 1), np.percentile(x[mask], 99)
        return np.clip((x - lo) / (hi - lo + 1e-9), 0, 1)

    Hn, En = _norm(H, roi), _norm(E, roi)
    mx, mn = rgb.max(2), rgb.min(2)
    tissue = roi & ~(mn > 0.80) & ~(mx < 0.20)        # not white slide, not black border
    n_tissue = int(tissue.sum())
    if n_tissue < 1000:
        raise RuntimeError("too little tissue inside ROI")

    thr = threshold_otsu(Hn[tissue])
    nuclei = tissue & (Hn > thr)
    stroma = tissue & (En > np.percentile(En[tissue], eosin_pct)) & (Hn < thr)

    res = {"box_color": color,
           "tissue_frac_of_ROI": round(n_tissue / float(roi.sum()), 3),
           "nuclei_fraction": round(float(nuclei.sum()) / n_tissue, 3),
           "stroma_fraction": round(float(stroma.sum()) / n_tissue, 3),
           "mean_haematoxylin": round(float(Hn[tissue].mean()), 3),
           "mean_eosin": round(float(En[tissue].mean()), 3)}
    res["H_over_E"] = round(res["mean_haematoxylin"] / (res["mean_eosin"] + 1e-9), 3)
    res["note"] = ("cellularity/stroma proxy, NOT calibrated purity; lower nuclei_fraction / "
                   "higher stroma_fraction indicates a stroma-dominated, lower-tumour-content section")
    if adata is not None:
        adata.uns["spatial_omics_he_purity"] = res
    return res
