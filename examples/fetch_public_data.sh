#!/usr/bin/env bash
# Fetch public slide-DNA-seq data (Zhao et al. Nature 2021/2022) for the
# "reproduce Nature 2022" flagship tutorial. BioProject: PRJNA768453.
#
# Requires: sra-tools (prefetch, fasterq-dump). Install via:
#   conda install -c bioconda sra-tools
#
# NOTE: these are large raw-read runs — run on a machine with disk + bandwidth,
# ideally the same HPC that runs the upstream pipeline (stage1-2). This script
# only DOWNLOADS; converting reads -> spot×bin matrix is the upstream pipeline,
# after which you call spatial_omics.io.from_pipeline(...).
set -euo pipefail

OUTDIR="${1:-public_data/slide_dna_seq}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Edit this list to the specific runs you want (find them on the SRA run selector
# for PRJNA768453). Example placeholders — REPLACE with real SRR accessions:
RUNS=(
  # SRRXXXXXXX   # mouse metastasis slide-DNA-seq array 1
  # SRRXXXXXXX   # human brain tumour slide-DNA-seq array
)

if [ ${#RUNS[@]} -eq 0 ]; then
  echo "No SRR accessions set. Open the SRA Run Selector for PRJNA768453,"
  echo "pick the slide-DNA-seq runs, and add their SRR IDs to the RUNS=() array."
  echo "  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA768453"
  exit 1
fi

for srr in "${RUNS[@]}"; do
  echo ">> prefetch $srr"
  prefetch "$srr"
  echo ">> fasterq-dump $srr"
  fasterq-dump --split-files --threads 8 "$srr" -O .
  pigz -p 8 "${srr}"_*.fastq || gzip "${srr}"_*.fastq
done

echo "Done. Next: run the upstream pipeline (stage1-2) on these FASTQs to produce"
echo "top50x50_fragments_with_cb.tsv.gz + matched_spot_barcodes.tsv, then:"
echo "  python -c \"import spatial_omics as so; ad=so.io.from_pipeline('frag.tsv.gz','matched.tsv'); ...\""
