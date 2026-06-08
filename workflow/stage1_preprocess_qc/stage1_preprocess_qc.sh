#!/usr/bin/env bash
set -uo pipefail
###############################################################################
# Stage 1 : spatial DNA-seq preprocessing + bulk-like QC
#
# NOTE on naming: an earlier copy of this script was called "atac_qc.sh".
#   The data is spatial DNA-seq (genomic DNA), NOT ATAC. The ATAC-style QC
#   metrics (MACS2 peaks / FRiP / TSS enrichment) are kept only as generic
#   sanity QC: for genomic DNA, TSS enrichment is expected to be ~FLAT (≈1),
#   the opposite of ATAC. The barcode/alignment/dedup/fragment logic is the
#   real product and is unchanged from the validated pipeline.
#
# Runs BOTH samples (520_520, 525_525), each in its own  <sample>_out/  dir.
# Critical Stage-3 inputs (spatial_prep/*.q30.*) are produced right after the
# q30 BAM, BEFORE the optional QC, and optional QC steps are guarded so a
# failure there never blocks the Stage-3 inputs.
###############################################################################

#########################  CONFIG — edit these  ###############################
SAMPLES=("520_520" "525_525")
DATA_DIR="."                                   # dir holding <sample>_1.fq.gz / _2.fq.gz
BOWTIE2_INDEX="/home/data/reference/index/bowtie/hg38/hg38"
GTF="/home/data/reference/C4/hg38.gtf"

THREADS=48
LINKER_SEQ="CAGTCATGTCATGAGCTA"
ME_SEQ="CTGTCTCTTATACACATCT"
SEQ_START=113
BC2_START=0;  BC2_END=8
BC1_START=26; BC1_END=34
BARCODE_ORDER="BC2BC1"
MIN_R1_LEN=20
MAPQ_Q30=30
MACS2_GENOME="hs"
MACS2_QVALUE="0.01"
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_sample () {
  local SAMPLE="$1"
  local R1="${DATA_DIR}/${SAMPLE}_1.fq.gz"
  local R2="${DATA_DIR}/${SAMPLE}_2.fq.gz"
  local OUT="${SAMPLE}_out"
  echo "==================== SAMPLE ${SAMPLE} ===================="
  [[ -f "$R1" ]] || { echo "[ERROR] missing $R1"; return 1; }
  [[ -f "$R2" ]] || { echo "[ERROR] missing $R2"; return 1; }
  mkdir -p "${OUT}/logs" "${OUT}/spatial_prep"

  # 1. linker-read filtering (keep pairs whose R1 contains the linker)
  echo "[${SAMPLE}] 1) bbduk linker filter"
  bbduk.sh in1="$R1" in2="$R2" \
    outm1="${OUT}/r1.linker.fq.gz" outm2="${OUT}/r2.linker.fq.gz" \
    k=18 restrictleft=85 hdist=3 skipr2=t literal="${LINKER_SEQ}" \
    threads="${THREADS}" stats="${OUT}/linker_filter.stats.txt" \
    > "${OUT}/logs/bbduk.out.log" 2> "${OUT}/logs/bbduk.err.log" || { echo "[ERROR] bbduk failed"; return 1; }

  # 2. trim ME / Tn5 read-through adapter
  echo "[${SAMPLE}] 2) cutadapt ME trim"
  cutadapt -j "${THREADS}" -a "${ME_SEQ}" -A "${ME_SEQ}" -m 20 \
    -o "${OUT}/r1.linker.MEtrim.fq.gz" -p "${OUT}/r2.linker.MEtrim.fq.gz" \
    "${OUT}/r1.linker.fq.gz" "${OUT}/r2.linker.fq.gz" \
    > "${OUT}/logs/cutadapt.log" 2>&1 || { echo "[ERROR] cutadapt failed"; return 1; }

  # 3. extract BC2+BC1 -> read name, truncate R1 to genomic insert
  echo "[${SAMPLE}] 3) barcode injection + R1 trim"
  python "${SCRIPT_DIR}/make_trimmed_r1_with_barcode_v2.py" \
    -i1 "${OUT}/r1.linker.MEtrim.fq.gz" -i2 "${OUT}/r2.linker.MEtrim.fq.gz" \
    -o1 "${OUT}/r1.trimmed.with_CB.fq.gz" -o2 "${OUT}/r2.trimmed.with_CB.fq.gz" \
    --seq_start "${SEQ_START}" \
    --bc2_start "${BC2_START}" --bc2_end "${BC2_END}" \
    --bc1_start "${BC1_START}" --bc1_end "${BC1_END}" \
    --barcode_order "${BARCODE_ORDER}" --min_r1_len "${MIN_R1_LEN}" \
    --barcode_whitelist_out "${OUT}/barcodes.whitelist.txt" \
    --barcode_counts_out   "${OUT}/barcode_counts.tsv" \
    --drop_stats_out       "${OUT}/barcode_drop_stats.tsv" \
    --preview_out          "${OUT}/barcode_preview.tsv" \
    > "${OUT}/logs/barcode_injection.log" 2>&1 || { echo "[ERROR] barcode injection failed"; return 1; }

  # 4. bowtie2 alignment + coord sort
  #    (tip: append  -X 2000  below if you want larger proper-pair fragments)
  echo "[${SAMPLE}] 4) bowtie2"
  bowtie2 -p "${THREADS}" -x "${BOWTIE2_INDEX}" \
    -1 "${OUT}/r1.trimmed.with_CB.fq.gz" -2 "${OUT}/r2.trimmed.with_CB.fq.gz" \
    2> "${OUT}/logs/bowtie2.log" \
  | samtools sort -@ "${THREADS}" -o "${OUT}/aln.sorted.bam" || { echo "[ERROR] bowtie2/sort failed"; return 1; }
  samtools index "${OUT}/aln.sorted.bam"

  # 5. name-sort -> fixmate -> pos-sort -> markdup
  echo "[${SAMPLE}] 5) fixmate + markdup"
  samtools sort -@ "${THREADS}" -n -o "${OUT}/aln.name.bam" "${OUT}/aln.sorted.bam"
  samtools fixmate -@ "${THREADS}" -m "${OUT}/aln.name.bam" "${OUT}/aln.fixmate.bam"
  samtools sort -@ "${THREADS}" -o "${OUT}/aln.fixmate.pos.bam" "${OUT}/aln.fixmate.bam"
  samtools markdup -@ "${THREADS}" "${OUT}/aln.fixmate.pos.bam" "${OUT}/aln.markdup.bam" || { echo "[ERROR] markdup failed"; return 1; }
  samtools index "${OUT}/aln.markdup.bam"
  samtools flagstat "${OUT}/aln.markdup.bam" > "${OUT}/flagstat.markdup.txt"

  # 6. q30 / properly-paired filter
  echo "[${SAMPLE}] 6) q30 filter"
  samtools view -@ "${THREADS}" -b -f 2 -F 1804 -q "${MAPQ_Q30}" \
    "${OUT}/aln.markdup.bam" > "${OUT}/aln.pp.q30.bam" || { echo "[ERROR] q30 filter failed"; return 1; }
  samtools index "${OUT}/aln.pp.q30.bam"
  samtools flagstat "${OUT}/aln.pp.q30.bam" > "${OUT}/flagstat.pp.q30.txt"

  ######  CRITICAL Stage-3 inputs — produced early, NOT guarded  ######
  echo "[${SAMPLE}] *) per-barcode fragment files (Stage-3 inputs)"
  # per-barcode high-quality fragment counts
  samtools view -f 67 "${OUT}/aln.pp.q30.bam" \
  | awk 'BEGIN{OFS="\t"}
    { if (match($1,/CB=[A-Z]+/)) { cb=substr($1,RSTART+3,RLENGTH-3); count[cb]++ } }
    END{ print "barcode","fragments"; for(b in count) print b,count[b] }' \
    > "${OUT}/spatial_prep/spot_fragments.q30.tsv" || { echo "[ERROR] spot_fragments failed"; return 1; }
  # per-fragment + barcode (chr,start,end,cb,mapq)
  samtools view -f 67 "${OUT}/aln.pp.q30.bam" \
  | awk 'BEGIN{OFS="\t"}
    { chr=$3; start=$4-1; len=$9; if(len<0)len=-len; end=start+len; mapq=$5;
      cb="NA"; if(match($1,/CB=[A-Z]+/)) cb=substr($1,RSTART+3,RLENGTH-3);
      print chr,start,end,cb,mapq }' \
  | gzip > "${OUT}/spatial_prep/fragments_with_cb.q30.tsv.gz" || { echo "[ERROR] fragments_with_cb failed"; return 1; }
  echo "[${SAMPLE}]    -> ${OUT}/spatial_prep/  (these feed Stage 2 then Stage 3)"

  ######  Optional bulk-like QC — guarded so failures don't block above  ######
  echo "[${SAMPLE}] 7-12) optional QC (guarded)"
  {
    # 7. fragment length distribution
    samtools view -@ "${THREADS}" -f 2 "${OUT}/aln.pp.q30.bam" \
    | awk '{t=$9; if(t<0)t=-t; if(t>0 && t<=1000) print t}' \
    | sort -n | uniq -c > "${OUT}/tlen.hist.q30.txt"

    # 8. peak calling (genomic DNA may yield few peaks; that is fine)
    macs2 callpeak -t "${OUT}/aln.pp.q30.bam" -f BAMPE -g "${MACS2_GENOME}" \
      -n bulk_q30 --outdir "${OUT}/macs2_bulk_q30" -q "${MACS2_QVALUE}" \
      > "${OUT}/logs/macs2.log" 2>&1

    # 9. FRiP
    if [[ -f "${OUT}/macs2_bulk_q30/bulk_q30_peaks.narrowPeak" ]]; then
      cut -f1-3 "${OUT}/macs2_bulk_q30/bulk_q30_peaks.narrowPeak" > "${OUT}/peaks.bed"
      TOTAL=$(samtools view -@ "${THREADS}" -c "${OUT}/aln.pp.q30.bam")
      INPEAK=$(bedtools intersect -u -abam "${OUT}/aln.pp.q30.bam" -b "${OUT}/peaks.bed" | samtools view -c -)
      { echo -e "TOTAL_READS_Q30\t${TOTAL}"; echo -e "IN_PEAK_READS_Q30\t${INPEAK}";
        awk -v a="${INPEAK}" -v b="${TOTAL}" 'BEGIN{print "FRiP_Q30\t" (b>0?a/b:0)}'; } > "${OUT}/frip.q30.txt"
    fi

    # 10. TSS bed from GTF
    if [[ -f "${GTF}" ]]; then
      awk 'BEGIN{FS=OFS="\t"} !/^#/ && ($3=="transcript"||$3=="mRNA"){
        chr=$1; if($7=="+"){t=$4-1;t2=$4} else if($7=="-"){t=$5-1;t2=$5} else next;
        if(t<0)t=0; print chr,t,t2 }' "${GTF}" | sort -k1,1 -k2,2n | uniq > "${OUT}/tss.bed"

      # 11. cutsite bed
      bedtools bamtobed -i "${OUT}/aln.pp.q30.bam" \
      | awk 'BEGIN{OFS="\t"}{if($6=="+"){s=$2;e=$2+1}else{s=$3-1;e=$3} print $1,s,e,$4,1,$6}' \
      > "${OUT}/cutsites.q30.bed"

      # 12. TSS enrichment (expected ~flat for genomic DNA)
      TSS_BED="${OUT}/tss.bed" CUT_BED="${OUT}/cutsites.q30.bed" OUT_DIR="${OUT}" python3 - <<'PY'
import os, bisect, numpy as np, pandas as pd
from collections import defaultdict
tss=pd.read_csv(os.environ["TSS_BED"],sep="\t",header=None,names=["chr","start","end"])
cs =pd.read_csv(os.environ["CUT_BED"],sep="\t",header=None,names=["chr","start","end","name","score","strand"])
cuts=defaultdict(list)
for c,s in zip(cs["chr"].values,cs["start"].values): cuts[c].append(int(s))
for c in cuts: cuts[c].sort()
bins=np.arange(-2000,2000,10); counts=np.zeros(len(bins),dtype=np.int64)
for c,pos in zip(tss["chr"].values,tss["start"].values):
    arr=cuts.get(c)
    if not arr: continue
    center=int(pos); i=bisect.bisect_left(arr,center-2000); j=bisect.bisect_right(arr,center+2000)
    for p in arr[i:j]:
        b=(p-center+2000)//10
        if 0<=b<len(counts): counts[b]+=1
df=pd.DataFrame({"offset":bins,"count":counts})
df.to_csv(os.path.join(os.environ["OUT_DIR"],"tss_profile_10bp.q30.tsv"),sep="\t",index=False)
peak=df[(df.offset>=-50)&(df.offset<=50)]["count"].mean()
flank=df[((df.offset>=-2000)&(df.offset<=-1500))|((df.offset>=1500)&(df.offset<=2000))]["count"].mean()
score=float(peak/flank) if flank>0 else float("nan")
open(os.path.join(os.environ["OUT_DIR"],"tss_enrichment.q30.txt"),"w").write(f"TSS_enrichment_q30\t{score}\n")
PY
    fi

    # 14. plot (optional helper)
    if [[ -f "${SCRIPT_DIR}/plot_bulk_atac_qc.py" ]]; then
      ( cd "${OUT}" && python3 "${SCRIPT_DIR}/plot_bulk_atac_qc.py" --logy )
    fi
  } || echo "[${SAMPLE}] [warn] optional QC step failed — Stage-3 inputs are already written, continuing."

  echo "[${SAMPLE}] DONE."
}

for S in "${SAMPLES[@]}"; do
  run_sample "$S" || echo "[ERROR] sample $S aborted (see logs)"
done
echo "[INFO] Stage 1 finished. Next: run stage2_spatial_map for each sample on its spatial_prep/ outputs."
