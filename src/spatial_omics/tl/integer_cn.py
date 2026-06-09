"""Integer copy number + purity/ploidy from segments.

Observed copy-number amplitude is compressed by tumour impurity: a spot that is fraction
``p`` tumour shows ``obs = p·CN_tumour + (1-p)·2`` (diploid normal contamination). So a real
single-copy gain (CN 3) at p=0.5 reads as 2.5, not 3. To assign INTEGER copy-number states
we must recover ``p`` (the scale). We grid-fit the purity that makes the segment means land
on an integer lattice — a BAF-free ASCAT/Ginkgo-style fit (we have no allele frequencies at
this depth). Degeneracy (p and p/2 both fit) is broken by parsimony: the highest purity among
near-optimal fits (least rescaling, smallest copy-number states).

Guard: if the genome is event-free, purity is unidentifiable from CNA (CN = 2N everywhere) →
flagged, not forced. If amplitudes are too compressed to snap to integers, ``resolvable`` is
False — treat integer CN as provisional and prefer an external purity (RNA / pathology).
Requires ``tl.segment`` first.
"""
from __future__ import annotations

import numpy as np


def integer_cn(adata, *, purity=None, baseline=2.0, max_cn=8, min_purity=0.1, n_grid=181,
               event_thresh=0.1):
    """Assign integer copy-number states to segments; fit (or take) tumour purity.

    Parameters
    ----------
    purity : float, optional
        Tumour fraction in (0, 1]. If None, it is fitted from the CNA segment lattice.
    baseline : float
        Observed value of a diploid (2N) segment (``copy_number`` anchors autosomes to 2).
    max_cn : int
        Clip integer states to [0, max_cn].
    min_purity, n_grid : float, int
        Purity grid search range/resolution when ``purity`` is None.
    event_thresh : float
        A segment is an "event" (constrains the fit) if |mean − baseline| exceeds this.

    Writes ``adata.var['integer_cn']`` (per-bin integer state; -1 unsegmented) and
    ``adata.uns['spatial_omics_integer_cn']`` (purity, ploidy, fit_residual, resolvable,
    per-segment table). Returns the segment table with ``integer_cn`` / ``cn_continuous``.
    """
    seg = adata.uns.get("spatial_omics_segments")
    if seg is None:
        raise KeyError("run tl.segment first")
    tbl = seg["table"].copy()
    means = tbl["mean"].values.astype(float)
    w = tbl["n_bins"].values.astype(float)
    dev = means - baseline
    has_events = bool(np.max(np.abs(dev)) > event_thresh)

    ev = np.abs(dev) > event_thresh             # event segments constrain the fit

    def lattice_cost(p):
        """weighted squared distance of segment CN to the integer lattice (NOT normalised,
        so the event segments — not the flat majority — drive the fit)."""
        cn = baseline + dev / p
        return float(np.sum(w * (cn - np.round(cn)) ** 2))

    if purity is not None:
        p = float(purity)
        fitted_p, source = p, "given"
    elif not has_events:
        p, fitted_p, source = 1.0, None, "unidentifiable (no events)"
    else:
        ps = np.linspace(min_purity, 1.0, n_grid)
        cost = np.array([lattice_cost(pp) for pp in ps])
        near = ps[cost <= cost.min() * 1.5 + 1e-6]   # all near-optimal fits (relative slack)
        p = float(near.max())                        # parsimony: highest purity / smallest states
        fitted_p, source = p, "fitted"

    cn_cont = baseline + dev / p
    cn_int = np.clip(np.round(cn_cont), 0, max_cn).astype(int)
    # fit quality = max distance of an EVENT segment from its integer (0 = perfect snap)
    max_int_dev = float(np.max(np.abs(cn_cont[ev] - np.round(cn_cont[ev])))) if ev.any() else 0.0

    seg_to_cn = dict(zip(tbl["segment"].values, cn_int))
    seg_id = adata.var["segment"].values
    bin_cn = np.array([seg_to_cn.get(s, -1) for s in seg_id], dtype=int)
    adata.var["integer_cn"] = bin_cn
    tbl["cn_continuous"] = np.round(cn_cont, 3)
    tbl["integer_cn"] = cn_int

    ploidy = float(np.average(cn_int, weights=w))
    resolvable = bool(has_events and max_int_dev < 0.15)
    if not has_events:
        verdict = "no CNA events: CN = 2N genome-wide; purity unidentifiable from copy number"
    elif resolvable:
        verdict = "integer CN resolved at purity %.2f (max integer deviation %.3f, %d states)" % (
            p, max_int_dev, len(set(cn_int)))
    else:
        verdict = ("amplitudes too compressed to snap to integers (low purity) — integer CN "
                   "provisional; prefer an external purity (RNA / pathology)")
    adata.uns["spatial_omics_integer_cn"] = {
        "purity": round(fitted_p, 3) if fitted_p is not None else None,
        "purity_source": source, "ploidy": round(ploidy, 2),
        "max_integer_deviation": round(max_int_dev, 4), "n_states": int(len(set(cn_int))),
        "resolvable": resolvable, "table": tbl, "verdict": verdict}
    return tbl
