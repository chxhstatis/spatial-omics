"""spatial_omics quickstart — runs the whole pipeline on synthetic data and saves figures.

    python examples/quickstart.py

This is also executed in CI as a smoke test, so the tutorial can never go stale.
"""
import os

import spatial_omics as sc

OUT = os.path.join(os.path.dirname(__file__), "_figures")
os.makedirs(OUT, exist_ok=True)

# 1. Get data — synthetic here; swap for sc.io.from_pipeline(frag, matched) on real data.
adata = sc.datasets.simulate(coverage=0.4, n_clones=2, seed=0)
print(adata)

# 2. Preprocess: QC bins, correct GC/mappability bias, densify + normalize.
sc.pp.bin_qc(adata)
sc.pp.correct_bias(adata)
sc.pp.normalize(adata)

# 3. Analyse: double-smooth PCs, call clones de-novo, derive copy number, test significance.
sc.tl.dual_smooth(adata)
sc.tl.call_clones(adata)
sc.tl.copy_number(adata)
sc.tl.permutation_significance(adata, regions=["chr8", "chr18"])

print("clones:", adata.uns["spatial_omics_clones"]["sizes"])

# 4. Plot.
sc.pl.spatial_clones(adata).savefig(f"{OUT}/clones.png", dpi=150, bbox_inches="tight")
sc.pl.spatial_copy_number(adata, chrom="chr8").savefig(f"{OUT}/cn_chr8.png", dpi=150, bbox_inches="tight")
sc.pl.significance(adata, region="chr8").savefig(f"{OUT}/sig_chr8.png", dpi=150, bbox_inches="tight")
sc.pl.clone_profiles(adata).savefig(f"{OUT}/clone_profiles.png", dpi=150, bbox_inches="tight")
print("figures ->", OUT)
