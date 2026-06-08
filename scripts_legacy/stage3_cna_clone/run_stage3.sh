#!/usr/bin/env bash
set -uo pipefail
###############################################################################
# Stage 3 runner (HPC-friendly, self-contained).
# Usage:
#   bash run_stage3.sh <stage2_dir> [outroot] [sample1 sample2 ...]
#   - <stage2_dir>: dir containing the Stage-2 outputs. Each sample may be in
#     <stage2_dir>/<sample>/  OR  <stage2_dir>/<sample>_xy_map/  (auto-detected).
#     Required files per sample: top50x50_fragments_with_cb.tsv.gz, matched_spot_barcodes.tsv
#   - [outroot]: output prefix (default: stage3_outputs)
#   - [samples]: default "520_520 525_525"
# Reference files (GC / mappability / blacklist) are bundled in ./ref — no downloads needed.
###############################################################################
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE2_DIR="${1:-.}"
OUTROOT="${2:-stage3_outputs}"
SAMPLES=("${@:3}")
if [ ${#SAMPLES[@]} -eq 0 ]; then SAMPLES=(520_520 525_525); fi

echo "[stage3] stage2_dir=${STAGE2_DIR}  outroot=${OUTROOT}  samples=${SAMPLES[*]}"

# (D) PRIMARY: bulk tumor CNA profile + spatial CNA-burden map (GC+mappability+blacklist)
python3 "${HERE}/stage3_bulk_cna.py" \
  --stage2_dir "${STAGE2_DIR}" --samples "${SAMPLES[@]}" \
  --ref_dir "${HERE}/ref" --outroot "${OUTROOT}_bulk" \
  && echo "[stage3] bulk CNA done -> ${OUTROOT}_bulk/" || echo "[stage3] bulk CNA FAILED"

# clone detection (GC + depth-confound corrected; may find no robust subclones at low purity)
python3 "${HERE}/stage3_spatial_clones.py" \
  --stage2_dir "${STAGE2_DIR}" --samples "${SAMPLES[@]}" \
  --gc_file "${HERE}/ref/bin_gc_1mb.tsv" --outroot "${OUTROOT}_clones" \
  --bin_size 1000000 --spatial_win 7 \
  && echo "[stage3] clone detection done -> ${OUTROOT}_clones/" || echo "[stage3] clone detection FAILED"

# PER-SPOT CNV: full genome-wide copy number for EACH of the 2500 spots (spatial smoothing)
python3 "${HERE}/stage3_per_spot_cnv.py" \
  --stage2_dir "${STAGE2_DIR}" --samples "${SAMPLES[@]}" \
  --ref_dir "${HERE}/ref" --outroot "${OUTROOT}_perspot" \
  --spatial_win 7 --cn_chroms chr8,chr17,chr18 \
  && echo "[stage3] per-spot CNV done -> ${OUTROOT}_perspot/" || echo "[stage3] per-spot CNV FAILED"

echo "[stage3] ALL DONE. Key outputs:"
echo "  ${OUTROOT}_bulk/<sample>/bulk_cna_profile.png  + spatial_cna_burden_map.png"
echo "  ${OUTROOT}_perspot/<sample>/per_spot_cnv_matrix.tsv.gz  (2500 spots x ~2646 bins) + spatial_cn_<chr>.png"
