"""Cohort comparison — cross-sample bulk CNA, with the sex (chrX) confound handled.

Compares several samples' pseudobulk copy-number profiles: which CNA events recur across
patients (drivers) vs are private, and which tumours are alike. Sample-similarity is
computed on AUTOSOMES ONLY — including chrX would cluster samples by sex, not tumour
biology (a real trap: two males look 0.9-correlated just from a shared chrX loss).
Ported from the stage3 pipeline's ``compare_samples_cna`` (original in git history /
the project archive); this AnnData-native version is now the single implementation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .._constants import autosomes, chrom_order

#: canonical PDAC interpretation loci (gene -> (chr, Mb))
_LOCI = {"CDKN2A(9p)": ("chr9", 21.0), "MYC(8q)": ("chr8", 127.0), "KRAS(12p)": ("chr12", 25.0),
         "TP53(17p)": ("chr17", 7.0), "SMAD4(18q)": ("chr18", 51.0), "GATA6(18q11)": ("chr18", 22.0)}


def _bulk_log2(adata, layer):
    """Per-bin pseudobulk log2, re-centred on this sample's autosomal median (so chrX
    reads ~-1 in males regardless of how ``copy_number`` anchored)."""
    cn = np.asarray(adata.layers[layer], dtype=float)
    tot = np.asarray(adata.layers["counts"]).sum(1)
    cn = cn[tot > 0]
    pb = np.nanmedian(cn, axis=0)
    auto = adata.var["chr"].isin(autosomes()).values
    base = np.nanmedian(pb[auto & np.isfinite(pb)])
    if not np.isfinite(base) or base <= 0:
        base = np.nanmedian(pb[np.isfinite(pb)])
    return np.log2(np.clip(pb / base, 1e-3, None))


def cohort_compare(adatas, *, layer="copy_number", thresh=0.2, male_chrx=-0.3):
    """Cross-sample bulk CNA comparison with the sex confound controlled.

    Parameters
    ----------
    adatas : dict[str, AnnData] | list[AnnData]
        Samples to compare. Each must have ``layers[layer]`` (run ``tl.copy_number``) and
        share the same bins (same genome + bin size).
    layer : str
        Per-spot copy-number layer to pseudobulk (default ``'copy_number'``).
    thresh : float
        |log2| gain/loss call threshold for the recurrence track.
    male_chrx : float
        chrX mean-log2 below which a sample is called male.

    Returns
    -------
    dict with:
        ``profiles`` — DataFrame (bins × samples) of bulk log2;
        ``correlation`` — autosomal sample×sample correlation (chrX EXCLUDED — the honest one);
        ``correlation_with_chrX`` — same including chrX (shows the sex confound);
        ``inferred_sex`` / ``chrX_mean_log2`` — per sample;
        ``recurrence`` — DataFrame(chr, start, end, n_gain, n_loss, mean_log2);
        ``recurrent_loci`` — canonical PDAC loci with per-sample log2 + cohort gain/loss counts.
    """
    if isinstance(adatas, dict):
        names, objs = list(adatas.keys()), list(adatas.values())
    else:
        objs = list(adatas)
        names = [a.uns.get("sample", f"sample{i}") for i, a in enumerate(objs)]

    ref = objs[0].var
    key = (ref["chr"].astype(str) + ":" + ref["start"].astype(str)).values
    prof = pd.DataFrame(index=key)
    prof["chr"] = ref["chr"].values
    prof["start"] = ref["start"].values
    prof["end"] = ref["end"].values if "end" in ref else ref["start"].values
    for name, a in zip(names, objs):
        akey = (a.var["chr"].astype(str) + ":" + a.var["start"].astype(str)).values
        prof[name] = pd.Series(_bulk_log2(a, layer), index=akey).reindex(key).values

    L = prof[names]
    auto = prof["chr"].isin(autosomes()).values
    is_x = (prof["chr"] == "chrX").values

    corr_auto = L[auto].corr()
    corr_all = L.corr()
    chrx_mean = {n: float(np.nanmean(L[n].values[is_x])) for n in names}
    sex = {n: ("male" if chrx_mean[n] < male_chrx else "female") for n in names}

    Lv = L.values
    n_gain = np.nansum(Lv > thresh, axis=1)
    n_loss = np.nansum(Lv < -thresh, axis=1)
    rec = prof[["chr", "start", "end"]].copy()
    rec["n_gain"] = n_gain.astype(int)
    rec["n_loss"] = n_loss.astype(int)
    rec["mean_log2"] = np.nanmean(Lv, axis=1)

    rows = []
    for locus, (chrom, mb) in _LOCI.items():
        pos = int(mb * 1e6)
        sel = (prof["chr"].values == chrom) & (prof["start"].values <= pos) & (prof["end"].values > pos)
        if not sel.any():
            continue
        vals = Lv[sel][0]
        row = {"locus": locus}
        row.update({n: round(float(v), 3) for n, v in zip(names, vals)})
        row["n_gain"] = int(np.nansum(vals > thresh))
        row["n_loss"] = int(np.nansum(vals < -thresh))
        rows.append(row)

    return {"samples": names, "profiles": prof,
            "correlation": corr_auto, "correlation_with_chrX": corr_all,
            "inferred_sex": sex, "chrX_mean_log2": {n: round(v, 3) for n, v in chrx_mean.items()},
            "recurrence": rec, "recurrent_loci": pd.DataFrame(rows)}
