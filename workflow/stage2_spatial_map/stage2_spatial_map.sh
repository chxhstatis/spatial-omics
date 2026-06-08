#!/usr/bin/env bash
set -uo pipefail
###############################################################################
# Stage 2 : DBiT spatial mapping (50x50 grid)
#
# For each sample, maps the 16 bp spot barcode (BC2+BC1) to its (x_id,y_id)
# coordinate on the 50x50 microfluidic grid, using the fixed chip design.
#   BC2 -> X axis (192-ref, select 50),  BC1 -> Y axis (192-ref, select 50),
#   <=1 mismatch per half, then crop to the selected 50x50 and remap to true
#   microfluidic coordinates.
#
# Input  (from Stage 1):  <sample>_out/spatial_prep/spot_fragments.q30.tsv
#                         <sample>_out/spatial_prep/fragments_with_cb.q30.tsv.gz
# Output (per sample):    <sample>_xy_map/  (mapped_spot_xy.tsv, top50x50_matrix.tsv,
#                         top50x50_fragments_with_cb.tsv.gz, heatmaps, ...)
###############################################################################

SAMPLES=("520_520" "525_525")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
X_REF="${SCRIPT_DIR}/X_8bp_192corrected.txt"
Y_REF="${SCRIPT_DIR}/Y_8bp_192corrected.txt"

for SAMPLE in "${SAMPLES[@]}"; do
  SP="${SAMPLE}_out/spatial_prep"
  OUT="${SAMPLE}_xy_map"
  echo "==================== Stage 2 : ${SAMPLE} ===================="
  if [[ ! -f "${SP}/spot_fragments.q30.tsv" ]]; then
    echo "[ERROR] missing ${SP}/spot_fragments.q30.tsv (run Stage 1 first)"; continue
  fi
  python "${SCRIPT_DIR}/fixed_xy50_keep_output.py" \
    --spot_fragments "${SP}/spot_fragments.q30.tsv" \
    --x_ref "${X_REF}" \
    --y_ref "${Y_REF}" \
    --fragments_with_cb "${SP}/fragments_with_cb.q30.tsv.gz" \
    --outdir "${OUT}" \
    --order XY \
    --max_mm_x 1 \
    --max_mm_y 1 \
  && echo "[${SAMPLE}] -> ${OUT}/  (mapped_spot_xy.tsv, top50x50_matrix.tsv, top50x50_fragments_with_cb.tsv.gz)" \
  || echo "[ERROR] Stage 2 failed for ${SAMPLE}"
done
echo "[INFO] Stage 2 finished. Copy each <sample>_xy_map/ back to the Mac for Stage 3."
