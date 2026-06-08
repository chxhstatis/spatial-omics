"""Run the full spatial_omics pipeline on REAL Stage-2 outputs and save results.

Usage:
    python examples/run_real_data.py \
        --fragments  /path/to/<sample>/top50x50_fragments_with_cb.tsv.gz \
        --matched    /path/to/<sample>/matched_spot_barcodes.tsv \
        --sample     520_520 \
        --outdir     results_520_520

Produces: <outdir>/<sample>.h5ad  +  clone / copy-number / significance figures.
Uses the packaged hg38 bias tracks by default (no download needed).
"""
import argparse
import os

import spatial_omics as so


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fragments", required=True, help="top50x50_fragments_with_cb.tsv.gz")
    ap.add_argument("--matched", required=True, help="matched_spot_barcodes.tsv")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--ref_dir", default=None, help="custom bias-track dir (else bundled)")
    ap.add_argument("--cn_chroms", default="chr8,chr18", help="comma list to paint as CN maps")
    ap.add_argument("--sig_regions", default="chr8,chr17,chr18", help="comma list for significance")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{args.sample}] loading Stage-2 outputs ...")
    adata = so.io.from_pipeline(args.fragments, args.matched, sample=args.sample)
    print(adata)

    so.pp.load_reference_tracks(adata, args.ref_dir)
    so.pp.bin_qc(adata)
    so.pp.correct_bias(adata)
    so.pp.normalize(adata)

    so.tl.dual_smooth(adata)
    so.tl.call_clones(adata)
    so.tl.copy_number(adata)
    so.tl.permutation_significance(adata, regions=args.sig_regions.split(","))

    print(f"[{args.sample}] clones:", adata.uns["spatial_omics_clones"]["sizes"])
    print(f"[{args.sample}] depth confound r(PC1):",
          round(adata.uns["spatial_omics_dual_smooth"]["depth_confound_r_pc1"], 3))

    h5 = os.path.join(args.outdir, f"{args.sample}.h5ad")
    adata.write(h5)
    so.pl.spatial_clones(adata).savefig(f"{args.outdir}/{args.sample}_clones.png", dpi=150, bbox_inches="tight")
    so.pl.clone_profiles(adata).savefig(f"{args.outdir}/{args.sample}_clone_profiles.png", dpi=150, bbox_inches="tight")
    for chrom in args.cn_chroms.split(","):
        so.pl.spatial_copy_number(adata, chrom=chrom).savefig(
            f"{args.outdir}/{args.sample}_cn_{chrom}.png", dpi=150, bbox_inches="tight")
    for region in args.sig_regions.split(","):
        so.pl.significance(adata, region=region).savefig(
            f"{args.outdir}/{args.sample}_sig_{region}.png", dpi=150, bbox_inches="tight")
    print(f"[{args.sample}] DONE -> {h5} (+ figures in {args.outdir}/)")


if __name__ == "__main__":
    main()
