"""spatial_omics — spatial copy-number & clone analysis for sparse spatial DNA-seq.

A scverse-style, AnnData-native toolkit for the data type slide-DNA-seq / DBiT
spatial DNA-seq produces: a spot x genomic-bin fragment matrix. Unlike inferCNV /
CopyKAT / Numbat (which infer CNV from RNA *expression*), spatial_omics works on the
*measured DNA reads* directly.

Typical flow::

    import spatial_omics as sc
    adata = sc.datasets.simulate(coverage=0.4, n_clones=2)   # or sc.io.from_pipeline(...)
    sc.pp.bin_qc(adata); sc.pp.correct_bias(adata); sc.pp.normalize(adata)
    sc.tl.dual_smooth(adata)
    sc.tl.call_clones(adata)
    sc.tl.copy_number(adata)
    sc.tl.permutation_significance(adata, regions=["chr8", "chr18"])
    sc.pl.spatial_clones(adata)

Submodules mirror scanpy: ``pp`` (preprocess), ``tl`` (tools), ``pl`` (plots),
plus ``io`` and ``datasets``.
"""
from __future__ import annotations

from . import datasets, io, pl, pp, tl
from .io import from_matrix, from_pipeline, read_h5ad

__version__ = "0.1.0"
__all__ = ["pp", "tl", "pl", "io", "datasets", "from_pipeline", "from_matrix", "read_h5ad"]
