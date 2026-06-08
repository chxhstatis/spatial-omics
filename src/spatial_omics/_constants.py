"""Genome constants shared across spatial_omics.

Currently bundles hg38 autosomes + chrX (chrY/chrM/contigs excluded; sex ignored),
matching the spatial_omics_pipeline reference design. New genomes (mm10, T2T) plug
in here by adding their chromosome-size dict and registering it in ``GENOMES``.
"""
from __future__ import annotations

HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895,
}

GENOMES = {"hg38": HG38}

#: chromosome plotting/order convention
def chrom_order(genome: str = "hg38") -> list[str]:
    return [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def autosomes(genome: str = "hg38") -> list[str]:
    return [f"chr{i}" for i in range(1, 23)]


#: PDAC interpretation anchors (gene -> (chr, Mb)) — used by plotting helpers.
PDAC_GENES = {
    "KRAS": ("chr12", 25.21), "CDKN2A": ("chr9", 21.97),
    "TP53": ("chr17", 7.67), "SMAD4": ("chr18", 51.03),
    "MYC": ("chr8", 127.74), "GATA6": ("chr18", 22.17),
}
