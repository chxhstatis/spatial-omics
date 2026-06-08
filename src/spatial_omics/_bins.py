"""Genomic bin construction (1 Mb tiling, matching stage3)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ._constants import GENOMES, chrom_order


def build_bins(bin_size: int = 1_000_000, genome: str = "hg38") -> pd.DataFrame:
    """Tile the genome into fixed-width bins.

    Returns a DataFrame with columns ``chr, start, end`` in genome order — this is
    the canonical column (``var``) order used everywhere in spatial_omics.
    """
    sizes = GENOMES[genome]
    rows = []
    for c in chrom_order(genome):
        n = (sizes[c] + bin_size - 1) // bin_size
        for b in range(n):
            rows.append((c, b * bin_size, min((b + 1) * bin_size, sizes[c])))
    bins = pd.DataFrame(rows, columns=["chr", "start", "end"])
    bins.index = bins["chr"] + ":" + bins["start"].astype(str)
    return bins


def bin_offsets(bin_size: int = 1_000_000, genome: str = "hg38") -> tuple[dict, int]:
    """Per-chromosome start offset into the flat bin array, and total #bins."""
    sizes = GENOMES[genome]
    off, gi = {}, 0
    for c in chrom_order(genome):
        off[c] = gi
        gi += (sizes[c] + bin_size - 1) // bin_size
    return off, gi
