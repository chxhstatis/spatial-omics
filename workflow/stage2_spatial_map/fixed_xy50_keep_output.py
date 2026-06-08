#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fixed_xy50_keep_output.py

Purpose
-------
Use pre-determined 50 X barcodes and 50 Y barcodes from 192-reference files,
then generate the same output files and plots as select_top50_xy_axis_keep_output.py,
without re-selecting top50 by fragment count.

Default selected reference indices, 1-based:
    X: 49-74,121-144  -> X1-X50
    Y: 49-70,73-78,123-144 -> Y1-Y50

Then remap X1-X50 / Y1-Y50 to Xnew1-Xnew50 / Ynew1-Ynew50 using the
microfluidic coordinate rule provided by the user. By default, x_id/y_id are
set to Xnew/Ynew, so downstream mapped_spot_xy.tsv uses true spatial order.
"""

import argparse
import gzip
import os
import re
from typing import Dict, Iterable, List, Tuple

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# Basic helpers
# ============================================================

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 999999
    return sum(x != y for x, y in zip(a, b))


def parse_index_ranges(spec: str) -> List[int]:
    """
    Parse strings like '49-74,121-144' into 1-based integer list.
    """
    out = []
    for part in str(spec).replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a, b = int(a), int(b)
            if b < a:
                raise ValueError(f"Bad range: {part}")
            out.extend(range(a, b + 1))
        else:
            out.append(int(part))
    if len(out) != len(set(out)):
        dup = [x for x in out if out.count(x) > 1]
        raise ValueError(f"Index ranges contain duplicated indices, e.g. {dup[:10]}")
    return out


def remap_coordinate_scalar(val: int) -> int:
    """
    User-provided mapping: loading order 1-50 -> true microfluidic order 1-50.
    """
    val = int(val)

    if 1 <= val <= 4:
        return val
    elif 47 <= val <= 50:
        return val - 42
    elif 5 <= val <= 8:
        return val + 4
    elif 43 <= val <= 46:
        return val - 30
    elif 9 <= val <= 12:
        return val + 8
    elif 39 <= val <= 42:
        return val - 18
    elif 13 <= val <= 16:
        return val + 12
    elif 35 <= val <= 38:
        return val - 6
    elif 17 <= val <= 19:
        return val + 16
    elif 32 <= val <= 34:
        return val + 4
    elif 20 <= val <= 22:
        return val + 19
    elif 29 <= val <= 31:
        return val + 13
    elif 23 <= val <= 25:
        return val + 22
    elif 26 <= val <= 28:
        return val + 22
    else:
        return val


def read_barcode_reference(path: str, axis: str) -> pd.DataFrame:
    """
    Read either:
      1) a plain 1-column 192-barcode file without header, or
      2) a TSV/CSV table containing X_barcode/Y_barcode.

    Returns columns:
        <axis>_ref_index, <axis>_barcode
    """
    axis = axis.upper()
    bc_col = f"{axis}_barcode"
    ref_col = f"{axis}_ref_index"

    with open(path, "rt") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        raise ValueError(f"Empty reference file: {path}")

    first_lower = lines[0].lower()
    looks_like_header = ("barcode" in first_lower) or ("coordinate" in first_lower) or ("ref_index" in first_lower)

    if looks_like_header:
        df = pd.read_csv(path, sep=None, engine="python")
        if bc_col not in df.columns:
            raise ValueError(f"{path} looks like a table but does not contain column {bc_col}. Current columns: {list(df.columns)}")

        out = pd.DataFrame()
        out[bc_col] = df[bc_col].astype(str).str.upper().str.strip()

        if ref_col in df.columns:
            out[ref_col] = df[ref_col].astype(int)
        else:
            # If no original reference index exists, use row order.
            out[ref_col] = range(1, len(out) + 1)

        out = out[[ref_col, bc_col]].copy()
    else:
        barcodes = []
        for ln in lines:
            # allow whitespace / tab / comma; use first field
            tok = re.split(r"[\t,\s]+", ln.strip())[0]
            if tok:
                barcodes.append(tok.upper())
        out = pd.DataFrame({
            ref_col: range(1, len(barcodes) + 1),
            bc_col: barcodes,
        })

    out[bc_col] = out[bc_col].astype(str).str.upper().str.strip()

    bad = out[out[bc_col].str.len() != 8]
    if bad.shape[0] > 0:
        raise ValueError(f"{path} contains non-8bp barcodes, e.g. {bad.head(5).to_dict(orient='records')}")

    if out[bc_col].duplicated().any():
        dup = out[out[bc_col].duplicated()][bc_col].tolist()
        raise ValueError(f"{path} contains duplicated barcodes, e.g. {dup[:10]}")

    return out


def build_fixed_axis(ref_df: pd.DataFrame, axis: str, selected_indices: List[int]) -> pd.DataFrame:
    """
    Build selected 50-axis table from 192 reference according to fixed original indices.

    Output columns:
        X_ref_index / Y_ref_index
        X_barcode / Y_barcode
        X_order / Y_order             loading order 1-50
        Xnew_order / Ynew_order       true microfluidic order 1-50
        X_coordinate / Y_coordinate   alias of Xnew/Ynew for compatibility
        x_id / y_id                   alias of Xnew/Ynew for output compatibility
    """
    axis = axis.upper()
    ref_col = f"{axis}_ref_index"
    bc_col = f"{axis}_barcode"
    order_col = f"{axis}_order"
    new_col = f"{axis}new_order"
    coord_col = f"{axis}_coordinate"

    ref_map = dict(zip(ref_df[ref_col], ref_df[bc_col]))
    missing = [idx for idx in selected_indices if idx not in ref_map]
    if missing:
        raise ValueError(f"{axis}: selected indices not found in reference file: {missing}")

    fixed = pd.DataFrame({
        ref_col: selected_indices,
        bc_col: [ref_map[i] for i in selected_indices],
        order_col: range(1, len(selected_indices) + 1),
    })

    if fixed.shape[0] != 50:
        raise ValueError(f"{axis}: expected 50 selected barcodes, got {fixed.shape[0]}")

    fixed[new_col] = fixed[order_col].map(remap_coordinate_scalar).astype(int)
    fixed[coord_col] = fixed[new_col].astype(int)

    # compatibility aliases used in old output files
    if axis == "X":
        fixed["x_id"] = fixed[new_col].astype(int)
        fixed["x_rank"] = fixed[new_col].astype(int)
    else:
        fixed["y_id"] = fixed[new_col].astype(int)
        fixed["y_rank"] = fixed[new_col].astype(int)

    if fixed[new_col].nunique() != 50 or set(fixed[new_col]) != set(range(1, 51)):
        raise ValueError(f"{axis}: remap_coordinate did not generate a complete 1-50 coordinate set")

    return fixed


def build_unique_best_matcher(refs: Iterable[str], max_mm: int = 1):
    """
    Build a fast unique-best matcher for 8bp barcode halves.

    For max_mm <= 1, this uses a precomputed exact/single-mismatch lookup, which is
    much faster than scanning all references for every observed barcode. For max_mm > 1,
    it falls back to the safe scan-based method.
    """
    refs = [str(r).upper().strip() for r in refs]
    ref_set = set(refs)

    if max_mm <= 0:
        def match_exact(query: str):
            query = str(query).upper().strip()
            if len(query) != 8:
                return "bad_length", None, None
            if query in ref_set:
                return "exact", query, 0
            return "unmapped", None, None
        return match_exact

    if max_mm == 1:
        alphabet = "ACGT"
        variant_to_ref: Dict[str, str] = {}
        ambiguous_variants = set()

        for ref in refs:
            for i, old_base in enumerate(ref):
                for base in alphabet:
                    if base == old_base:
                        continue
                    var = ref[:i] + base + ref[i + 1:]
                    if var in variant_to_ref and variant_to_ref[var] != ref:
                        ambiguous_variants.add(var)
                    else:
                        variant_to_ref[var] = ref

        def match_mm1(query: str):
            query = str(query).upper().strip()
            if len(query) != 8:
                return "bad_length", None, None
            if query in ref_set:
                return "exact", query, 0

            # Fast path for normal A/C/G/T barcodes.
            if set(query).issubset(set(alphabet)):
                if query in ambiguous_variants:
                    return "ambiguous", None, 1
                hit = variant_to_ref.get(query)
                if hit is not None:
                    return "corrected", hit, 1
                return "unmapped", None, None

            # Safe path for rare observed barcodes containing N or other characters.
            best_hits = [r for r in refs if hamming(query, r) <= 1]
            if len(best_hits) == 0:
                return "unmapped", None, None
            if len(best_hits) > 1:
                return "ambiguous", None, 1
            return "corrected", best_hits[0], 1

        return match_mm1

    # Rare fallback if someone explicitly uses max_mm > 1.
    def match_scan(query: str):
        query = str(query).upper().strip()
        if len(query) != 8:
            return "bad_length", None, None
        if query in ref_set:
            return "exact", query, 0
        best_mm = None
        best_hits = []
        for r in refs:
            mm = hamming(query, r)
            if mm <= max_mm:
                if best_mm is None or mm < best_mm:
                    best_mm = mm
                    best_hits = [r]
                elif mm == best_mm:
                    best_hits.append(r)
        if best_mm is None:
            return "unmapped", None, None
        if len(best_hits) > 1:
            return "ambiguous", None, best_mm
        return "corrected", best_hits[0], best_mm

    return match_scan


# ============================================================
# Plot helpers: keep old file names and general layout
# ============================================================

def plot_marginal(full_df, sel_df, coord_col, out_png, out_pdf, title):
    plt.figure(figsize=(8, 4))
    plt.bar(full_df[coord_col], full_df["fragments"])
    for v in sel_df[coord_col]:
        plt.axvline(v, linestyle="--", alpha=0.25)
    plt.xlabel(coord_col.upper())
    plt.ylabel("Total fragments")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.savefig(out_pdf)
    plt.close()


def plot_heatmap(df, x_order, y_order, out_png, out_pdf, title):
    x_to_i = {x: i for i, x in enumerate(x_order)}
    y_to_j = {y: j for j, y in enumerate(y_order)}
    mat = [[0 for _ in y_order] for _ in x_order]

    for _, row in df.iterrows():
        x = int(row["x_rank"])
        y = int(row["y_rank"])
        if x in x_to_i and y in y_to_j:
            mat[x_to_i[x]][y_to_j[y]] = row["fragments"]

    plt.figure(figsize=(10, 8))
    plt.imshow(mat, aspect="auto", origin="lower")
    plt.colorbar(label="Fragments")
    plt.xticks(range(len(y_order)), y_order, rotation=90, fontsize=6)
    plt.yticks(range(len(x_order)), x_order, fontsize=6)
    plt.xlabel("Ynew coordinate")
    plt.ylabel("Xnew coordinate")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.savefig(out_pdf)
    plt.close()


def plot_heatmap_by_coord(df, x_coords, y_coords, out_png, out_pdf, title):
    x_to_i = {x: i for i, x in enumerate(x_coords)}
    y_to_j = {y: j for j, y in enumerate(y_coords)}
    mat = [[0 for _ in y_coords] for _ in x_coords]

    for _, row in df.iterrows():
        x = int(row["x_id"])
        y = int(row["y_id"])
        if x in x_to_i and y in y_to_j:
            mat[x_to_i[x]][y_to_j[y]] = row["fragments"]

    plt.figure(figsize=(10, 8))
    plt.imshow(mat, aspect="auto", origin="lower")
    plt.colorbar(label="Fragments")
    plt.xticks(range(len(y_coords)), y_coords, rotation=90, fontsize=6)
    plt.yticks(range(len(x_coords)), x_coords, fontsize=6)
    plt.xlabel("Ynew coordinate")
    plt.ylabel("Xnew coordinate")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.savefig(out_pdf)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Match 16bp spot barcodes to fixed selected X/Y barcodes, "
            "use known Xnew/Ynew coordinates, and generate the same output files "
            "as select_top50_xy_axis_keep_output.py without selecting top50 again."
        )
    )

    ap.add_argument("--spot_fragments", required=True,
                    help="spatial_prep/spot_fragments.q30.tsv; must contain barcode and fragments/count columns")
    ap.add_argument("--x_ref", "--x_axis", dest="x_ref", required=True,
                    help="X_8bp_192corrected.txt, or a table containing X_barcode")
    ap.add_argument("--y_ref", "--y_axis", dest="y_ref", required=True,
                    help="Y_8bp_192corrected.txt, or a table containing Y_barcode")
    ap.add_argument("--fragments_with_cb", default="",
                    help="Optional: spatial_prep/fragments_with_cb.q30.tsv.gz")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--order", choices=["XY", "YX"], default="XY",
                    help="XY: barcode[:8]=X, barcode[8:]=Y; YX: barcode[:8]=Y, barcode[8:]=X")
    ap.add_argument("--max_mm_x", type=int, default=1)
    ap.add_argument("--max_mm_y", type=int, default=1)
    ap.add_argument("--x_ref_indices", default="49-74,121-144",
                    help="1-based X indices selected from the X reference file")
    ap.add_argument("--y_ref_indices", default="49-70,73-78,123-144",
                    help="1-based Y indices selected from the Y reference file")

    args = ap.parse_args()
    ensure_dir(args.outdir)

    # ============================================================
    # 1. Read 192 references and build fixed 50-axis tables
    # ============================================================
    x_ref = read_barcode_reference(args.x_ref, "X")
    y_ref = read_barcode_reference(args.y_ref, "Y")

    x_fixed = build_fixed_axis(x_ref, "X", parse_index_ranges(args.x_ref_indices))
    y_fixed = build_fixed_axis(y_ref, "Y", parse_index_ranges(args.y_ref_indices))

    x_fixed.to_csv(os.path.join(args.outdir, "fixed_X_ref_index_to_X1_50_and_Xnew1_50.tsv"), sep="\t", index=False)
    y_fixed.to_csv(os.path.join(args.outdir, "fixed_Y_ref_index_to_Y1_50_and_Ynew1_50.tsv"), sep="\t", index=False)

    # Match against ALL 192 references first, then crop to the known selected 50 x 50.
    # This avoids incorrectly correcting a sequencing-error barcode to a selected barcode
    # when its best match is actually a non-selected reference barcode.
    x_refs = x_ref["X_barcode"].tolist()
    y_refs = y_ref["Y_barcode"].tolist()
    x_matcher = build_unique_best_matcher(x_refs, max_mm=args.max_mm_x)
    y_matcher = build_unique_best_matcher(y_refs, max_mm=args.max_mm_y)

    # Full 192 reference annotations and selected-50 coordinate annotations.
    x_full_anno = x_ref.set_index("X_barcode").to_dict(orient="index")
    y_full_anno = y_ref.set_index("Y_barcode").to_dict(orient="index")
    x_selected_anno = x_fixed.set_index("X_barcode").to_dict(orient="index")
    y_selected_anno = y_fixed.set_index("Y_barcode").to_dict(orient="index")

    # ============================================================
    # 2. Read spot_fragments
    # ============================================================
    spot = pd.read_csv(args.spot_fragments, sep="\t")
    if "barcode" not in spot.columns:
        raise ValueError("spot_fragments file must contain a 'barcode' column")
    if "fragments" not in spot.columns:
        if "count" in spot.columns:
            spot = spot.rename(columns={"count": "fragments"})
        else:
            raise ValueError("spot_fragments file must contain 'fragments' or 'count' column")

    spot["barcode"] = spot["barcode"].astype(str).str.upper().str.strip()
    spot["fragments"] = spot["fragments"].astype(int)

    # ============================================================
    # 3. Split 16bp barcode and match only to fixed selected X/Y
    # ============================================================
    records = []

    for bc, fr in spot[["barcode", "fragments"]].itertuples(index=False, name=None):
        bc = str(bc).upper().strip()
        fr = int(fr)

        rec = {
            "barcode_obs": bc,
            "fragments": fr,
            "len": len(bc),
            "x_half_obs": None,
            "y_half_obs": None,
            "x_status": None,
            "y_status": None,
            "x_barcode": None,
            "y_barcode": None,
            "mm_x": None,
            "mm_y": None,
            "full_status": None,
            "barcode_mapped": None,
            "X_ref_index": None,
            "Y_ref_index": None,
            "X_order": None,
            "Y_order": None,
            "Xnew_order": None,
            "Ynew_order": None,
            # Compatibility with old downstream code: x_id/y_id are true Xnew/Ynew coordinates
            "x_id": None,
            "y_id": None,
            "x_rank": None,
            "y_rank": None,
        }

        if len(bc) != 16:
            rec["full_status"] = "bad_length"
            records.append(rec)
            continue

        if args.order == "XY":
            x_half = bc[:8]
            y_half = bc[8:]
        else:
            y_half = bc[:8]
            x_half = bc[8:]

        rec["x_half_obs"] = x_half
        rec["y_half_obs"] = y_half

        x_status, x_match, mm_x = x_matcher(x_half)
        y_status, y_match, mm_y = y_matcher(y_half)

        rec["x_status"] = x_status
        rec["y_status"] = y_status
        rec["x_barcode"] = x_match
        rec["y_barcode"] = y_match
        rec["mm_x"] = mm_x
        rec["mm_y"] = mm_y

        if x_status in {"exact", "corrected"} and y_status in {"exact", "corrected"}:
            rec["full_status"] = "mapped"
            rec["barcode_mapped"] = f"{x_match}{y_match}" if args.order == "XY" else f"{y_match}{x_match}"

            # Always record original 192-reference indices.
            rec["X_ref_index"] = int(x_full_anno[x_match]["X_ref_index"])
            rec["Y_ref_index"] = int(y_full_anno[y_match]["Y_ref_index"])

            # Only selected 50 x 50 barcodes receive X_order/Y_order and Xnew/Ynew coordinates.
            xa = x_selected_anno.get(x_match)
            ya = y_selected_anno.get(y_match)
            if xa is not None:
                rec["X_order"] = int(xa["X_order"])
                rec["Xnew_order"] = int(xa["Xnew_order"])
                rec["x_id"] = int(xa["Xnew_order"])
                rec["x_rank"] = int(xa["Xnew_order"])
            if ya is not None:
                rec["Y_order"] = int(ya["Y_order"])
                rec["Ynew_order"] = int(ya["Ynew_order"])
                rec["y_id"] = int(ya["Ynew_order"])
                rec["y_rank"] = int(ya["Ynew_order"])
        elif x_status == "ambiguous" or y_status == "ambiguous":
            rec["full_status"] = "ambiguous"
        elif x_status == "bad_length" or y_status == "bad_length":
            rec["full_status"] = "bad_length"
        else:
            rec["full_status"] = "unmapped"

        records.append(rec)

    match_df = pd.DataFrame(records)
    match_df.to_csv(os.path.join(args.outdir, "matched_spot_barcodes.tsv"), sep="\t", index=False)

    mapped = match_df[match_df["full_status"] == "mapped"].copy()

    # ============================================================
    # 4. X/Y marginal fragments over fixed 50 axes
    # ============================================================
    # Marginal counts are calculated for all matched 192-reference barcodes.
    # Selected 50 annotations are kept when available.
    x_marg = (
        mapped.groupby(["x_barcode", "X_ref_index"], as_index=False)["fragments"]
        .sum()
        .sort_values("X_ref_index")
    )
    x_marg = x_marg.merge(
        x_fixed[["X_barcode", "X_order", "Xnew_order"]].rename(columns={"X_barcode": "x_barcode"}),
        on="x_barcode", how="left"
    )
    x_marg["x_id"] = x_marg["Xnew_order"]
    x_marg["x_rank"] = x_marg["Xnew_order"]

    y_marg = (
        mapped.groupby(["y_barcode", "Y_ref_index"], as_index=False)["fragments"]
        .sum()
        .sort_values("Y_ref_index")
    )
    y_marg = y_marg.merge(
        y_fixed[["Y_barcode", "Y_order", "Ynew_order"]].rename(columns={"Y_barcode": "y_barcode"}),
        on="y_barcode", how="left"
    )
    y_marg["y_id"] = y_marg["Ynew_order"]
    y_marg["y_rank"] = y_marg["Ynew_order"]

    # Add missing fixed axes with zero fragments, so outputs/plots always have 50 rows.
    top_x = x_fixed[["X_ref_index", "X_barcode", "X_order", "Xnew_order"]].copy()
    top_x = top_x.rename(columns={"X_barcode": "x_barcode", "Xnew_order": "x_rank"})
    top_x["x_id"] = top_x["x_rank"].astype(int)
    top_x = top_x.merge(
        x_marg[["x_barcode", "fragments"]], on="x_barcode", how="left"
    )
    top_x["fragments"] = top_x["fragments"].fillna(0).astype(int)
    top_x = top_x.sort_values("x_id").reset_index(drop=True)

    top_y = y_fixed[["Y_ref_index", "Y_barcode", "Y_order", "Ynew_order"]].copy()
    top_y = top_y.rename(columns={"Y_barcode": "y_barcode", "Ynew_order": "y_rank"})
    top_y["y_id"] = top_y["y_rank"].astype(int)
    top_y = top_y.merge(
        y_marg[["y_barcode", "fragments"]], on="y_barcode", how="left"
    )
    top_y["fragments"] = top_y["fragments"].fillna(0).astype(int)
    top_y = top_y.sort_values("y_id").reset_index(drop=True)

    # For compatibility, output these filenames exactly.
    # all_matched_* contains all 192 references that received mapped fragments; selected rows have Xnew/Ynew annotations.
    x_marg_out = x_marg[["x_barcode", "X_ref_index", "X_order", "Xnew_order", "x_id", "fragments"]].copy()
    y_marg_out = y_marg[["y_barcode", "Y_ref_index", "Y_order", "Ynew_order", "y_id", "fragments"]].copy()

    x_marg_out.to_csv(os.path.join(args.outdir, "all_matched_x_fragments.tsv"), sep="\t", index=False)
    y_marg_out.to_csv(os.path.join(args.outdir, "all_matched_y_fragments.tsv"), sep="\t", index=False)
    top_x.to_csv(os.path.join(args.outdir, "top50_x.tsv"), sep="\t", index=False)
    top_y.to_csv(os.path.join(args.outdir, "top50_y.tsv"), sep="\t", index=False)

    # ============================================================
    # 5. Fixed 50 X × fixed 50 Y spot table
    # ============================================================
    cropped = mapped[
        mapped["x_barcode"].isin(set(x_fixed["X_barcode"])) &
        mapped["y_barcode"].isin(set(y_fixed["Y_barcode"]))
    ].copy()

    cropped_spot = (
        cropped.groupby(
            [
                "barcode_mapped",
                "x_barcode", "y_barcode",
                "X_ref_index", "Y_ref_index",
                "X_order", "Y_order",
                "Xnew_order", "Ynew_order",
                "x_id", "y_id", "x_rank", "y_rank",
            ],
            as_index=False,
        )["fragments"]
        .sum()
        .sort_values(["x_id", "y_id"])
    )

    cropped_spot.to_csv(os.path.join(args.outdir, "top50x50_spot_fragments.tsv"), sep="\t", index=False)

    mapped_spot_xy = (
        cropped_spot[["barcode_mapped", "x_id", "y_id"]]
        .drop_duplicates()
        .sort_values(["x_id", "y_id"])
    )
    mapped_spot_xy.to_csv(os.path.join(args.outdir, "mapped_spot_xy.tsv"), sep="\t", index=False)

    mat = (
        cropped_spot.pivot_table(index="x_id", columns="y_id", values="fragments", aggfunc="sum", fill_value=0)
    )
    mat = mat.reindex(index=list(range(1, 51)), columns=list(range(1, 51)), fill_value=0)
    mat.to_csv(os.path.join(args.outdir, "top50x50_matrix.tsv"), sep="\t")

    # Additional matrices useful for checking loading order vs true order.
    mat_loading = (
        cropped_spot.pivot_table(index="X_order", columns="Y_order", values="fragments", aggfunc="sum", fill_value=0)
    )
    mat_loading = mat_loading.reindex(index=list(range(1, 51)), columns=list(range(1, 51)), fill_value=0)
    mat_loading.to_csv(os.path.join(args.outdir, "top50x50_matrix_loading_X1Y1_order.tsv"), sep="\t")

    # ============================================================
    # 6. Summary
    # ============================================================
    summary = {
        "total_spot_barcodes_input": int(spot.shape[0]),
        "total_fragments_input": int(spot["fragments"].sum()),
        "fixed_x_n": int(top_x.shape[0]),
        "fixed_y_n": int(top_y.shape[0]),
        "mapped_spot_barcodes_all_refs": int(mapped.shape[0]),
        "mapped_fragments_all_refs": int(mapped["fragments"].sum()) if mapped.shape[0] else 0,
        "cropped_spot_barcodes": int(cropped_spot.shape[0]),
        "cropped_fragments": int(cropped_spot["fragments"].sum()) if cropped_spot.shape[0] else 0,
        "mean_fragments_per_cropped_spot": float(cropped_spot["fragments"].mean()) if cropped_spot.shape[0] else 0,
        "median_fragments_per_cropped_spot": float(cropped_spot["fragments"].median()) if cropped_spot.shape[0] else 0,
        "order": args.order,
        "max_mm_x": int(args.max_mm_x),
        "max_mm_y": int(args.max_mm_y),
        "x_ref_indices": args.x_ref_indices,
        "y_ref_indices": args.y_ref_indices,
    }

    with open(os.path.join(args.outdir, "summary.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        for k, v in summary.items():
            f.write(f"{k}\t{v}\n")

    # ============================================================
    # 7. Plots: same output names
    # ============================================================
    if top_x.shape[0] > 0:
        plot_marginal(
            top_x.rename(columns={"x_id": "coord"}),
            top_x.rename(columns={"x_id": "coord"}),
            "coord",
            os.path.join(args.outdir, "01_x_marginal_top50.png"),
            os.path.join(args.outdir, "01_x_marginal_top50.pdf"),
            "Fixed X barcodes: total fragments in Xnew order",
        )

    if top_y.shape[0] > 0:
        plot_marginal(
            top_y.rename(columns={"y_id": "coord"}),
            top_y.rename(columns={"y_id": "coord"}),
            "coord",
            os.path.join(args.outdir, "02_y_marginal_top50.png"),
            os.path.join(args.outdir, "02_y_marginal_top50.pdf"),
            "Fixed Y barcodes: total fragments in Ynew order",
        )

    if cropped_spot.shape[0] > 0:
        plot_heatmap(
            cropped_spot,
            list(range(1, 51)),
            list(range(1, 51)),
            os.path.join(args.outdir, "03_top50x50_heatmap_by_rank.png"),
            os.path.join(args.outdir, "03_top50x50_heatmap_by_rank.pdf"),
            "Fixed 50 X × fixed 50 Y heatmap (Xnew/Ynew order)",
        )

        plot_heatmap_by_coord(
            cropped_spot,
            list(range(1, 51)),
            list(range(1, 51)),
            os.path.join(args.outdir, "04_top50x50_heatmap_by_barcode_id.png"),
            os.path.join(args.outdir, "04_top50x50_heatmap_by_barcode_id.pdf"),
            "Fixed 50 X × fixed 50 Y heatmap (true coordinate order)",
        )

    # ============================================================
    # 8. Optional: filter fragment-level file
    # ============================================================
    if args.fragments_with_cb:
        selected_mapped_barcodes = set(cropped_spot["barcode_mapped"].tolist())
        mapped_match = match_df[
            (match_df["full_status"] == "mapped") &
            (match_df["barcode_mapped"].isin(selected_mapped_barcodes))
        ].copy()
        selected_obs_barcodes = set(mapped_match["barcode_obs"].astype(str).str.upper())

        out_frag_1 = os.path.join(args.outdir, "mapped_fragments_with_cb.tsv.gz")
        out_frag_2 = os.path.join(args.outdir, "top50x50_fragments_with_cb.tsv.gz")

        kept_n = 0
        total_frag_lines = 0
        with gzip.open(args.fragments_with_cb, "rt") as fin, \
             gzip.open(out_frag_1, "wt") as fout1, \
             gzip.open(out_frag_2, "wt") as fout2:
            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                total_frag_lines += 1
                cb = parts[3].upper()
                if cb in selected_obs_barcodes:
                    fout1.write(line + "\n")
                    fout2.write(line + "\n")
                    kept_n += 1

        with open(os.path.join(args.outdir, "fragment_filter_summary.tsv"), "w") as f:
            f.write("metric\tvalue\n")
            f.write(f"input_fragment_lines\t{total_frag_lines}\n")
            f.write(f"kept_fragments_lines\t{kept_n}\n")
            f.write(f"selected_obs_barcodes\t{len(selected_obs_barcodes)}\n")
            f.write(f"selected_mapped_barcodes\t{len(selected_mapped_barcodes)}\n")

    print("[INFO] Done.")
    print("[INFO] Matched against all references, then used fixed selected axes; no top50 re-selection was performed.")
    print("[INFO] Outputs:")
    print("  fixed_X_ref_index_to_X1_50_and_Xnew1_50.tsv")
    print("  fixed_Y_ref_index_to_Y1_50_and_Ynew1_50.tsv")
    print("  matched_spot_barcodes.tsv")
    print("  mapped_spot_xy.tsv")
    print("  top50_x.tsv")
    print("  top50_y.tsv")
    print("  all_matched_x_fragments.tsv")
    print("  all_matched_y_fragments.tsv")
    print("  top50x50_spot_fragments.tsv")
    print("  top50x50_matrix.tsv")
    print("  top50x50_matrix_loading_X1Y1_order.tsv")
    print("  summary.tsv")
    print("  01_x_marginal_top50.png/pdf")
    print("  02_y_marginal_top50.png/pdf")
    print("  03_top50x50_heatmap_by_rank.png/pdf")
    print("  04_top50x50_heatmap_by_barcode_id.png/pdf")
    if args.fragments_with_cb:
        print("  mapped_fragments_with_cb.tsv.gz")
        print("  top50x50_fragments_with_cb.tsv.gz")
        print("  fragment_filter_summary.tsv")


if __name__ == "__main__":
    main()
