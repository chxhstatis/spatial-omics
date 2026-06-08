#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_trimmed_r1_with_barcode_v2.py  (reconstructed)

Stage-1 step 3 of the spatial DNA-seq (DBiT 50x50) pipeline.

Input : linker-filtered, ME-trimmed paired FASTQ (R1, R2).
For every read pair it:
  1. extracts BC2 = R1[bc2_start:bc2_end]  and  BC1 = R1[bc1_start:bc1_end]
     (fixed positions; default BC2=R1[0:8], BC1=R1[26:34]),
  2. composes the 16 bp spatial barcode (BC2+BC1 by default; --barcode_order),
  3. injects the barcode into BOTH read names as 'CB=<barcode>' so the downstream
     awk can recover it with  match($1, /CB=[A-Z]+/)  ,
  4. truncates R1 to the genomic insert starting at --seq_start (default 113) and
     keeps R2 unchanged,
  5. drops a pair whose remaining R1 genomic length < --min_r1_len (default 20),
  6. writes a new paired FASTQ plus barcode whitelist / counts / drop-stats / preview.

Read-name format produced (so QNAME in the BAM carries the barcode):
    @<original_base>_CB=<BC2BC1>/1     and     @<original_base>_CB=<BC2BC1>/2
After alignment the aligner strips the trailing /1 /2, leaving the shared QNAME
'<original_base>_CB=<BC2BC1>', from which step-13 awk extracts the 16 bp CB.

Faithful default reproduces the original CLI:
  --seq_start 113 --bc2_start 0 --bc2_end 8 --bc1_start 26 --bc1_end 34
  --barcode_order BC2BC1 --min_r1_len 20

Optional robustness:
  --anchor_linker (+ --linker_seq) re-anchors BC1/BC2 (and the genomic start) to the
  ACTUAL linker position found per read, which fixes the ~9% of reads whose linker is
  shifted by an indel. OFF by default to match the validated fixed-position behaviour.
"""

import argparse
import gzip
import sys
from collections import Counter


def open_in(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def open_out_gz(path):
    # compresslevel=4 = good speed/size tradeoff for large FASTQ
    return gzip.open(path, "wt", compresslevel=4)


def name_base_and_mate(header):
    """'@NAME/1' or '@NAME comment' -> ('NAME', '/1'|'')."""
    h = header[1:] if header.startswith("@") else header
    h = h.rstrip("\n")
    # drop any Illumina-style comment after whitespace
    for sep in (" ", "\t"):
        if sep in h:
            h = h.split(sep, 1)[0]
    mate = ""
    if h.endswith("/1") or h.endswith("/2"):
        mate, h = h[-2:], h[:-2]
    return h, mate


def find_linker(seq, linker, expected_start, max_shift):
    """Return the 0-based start of the linker near expected_start (<=1 mismatch), or -1."""
    Llen = len(linker)
    lo = max(0, expected_start - max_shift)
    hi = expected_start + max_shift
    idx = seq.find(linker, lo, hi + Llen)
    if idx != -1:
        return idx
    best_mm, best_pos = 2, -1
    for cand in range(lo, hi + 1):
        seg = seq[cand:cand + Llen]
        if len(seg) < Llen:
            break
        mm = sum(1 for a, b in zip(seg, linker) if a != b)
        if mm <= 1 and mm < best_mm:
            best_mm, best_pos = mm, cand
    return best_pos


def main():
    ap = argparse.ArgumentParser(description="Inject spatial barcode into read names and trim R1 to genomic insert.")
    ap.add_argument("-i1", required=True, help="input R1 fastq(.gz)")
    ap.add_argument("-i2", required=True, help="input R2 fastq(.gz)")
    ap.add_argument("-o1", required=True, help="output R1 fastq.gz (genomic + CB in name)")
    ap.add_argument("-o2", required=True, help="output R2 fastq.gz (genomic + CB in name)")
    ap.add_argument("--seq_start", type=int, default=113, help="0-based start of genomic insert in R1")
    ap.add_argument("--bc2_start", type=int, default=0)
    ap.add_argument("--bc2_end", type=int, default=8)
    ap.add_argument("--bc1_start", type=int, default=26)
    ap.add_argument("--bc1_end", type=int, default=34)
    ap.add_argument("--barcode_order", choices=["BC2BC1", "BC1BC2"], default="BC2BC1")
    ap.add_argument("--min_r1_len", type=int, default=20, help="drop pair if trimmed R1 genomic < this")
    ap.add_argument("--barcode_whitelist_out", default="")
    ap.add_argument("--barcode_counts_out", default="")
    ap.add_argument("--drop_stats_out", default="")
    ap.add_argument("--preview_out", default="")
    ap.add_argument("--preview_n", type=int, default=100)
    # optional robust mode
    ap.add_argument("--anchor_linker", action="store_true",
                    help="re-anchor barcodes/genomic to the found linker position per read")
    ap.add_argument("--linker_seq", default="CAGTCATGTCATGAGCTA")
    ap.add_argument("--linker_max_shift", type=int, default=3)
    args = ap.parse_args()

    bc2_len = args.bc2_end - args.bc2_start
    bc1_len = args.bc1_end - args.bc1_start
    linker = args.linker_seq.upper()
    Llen = len(linker)

    bc_counts = Counter()
    drop = Counter()
    n_total = n_kept = 0
    preview = []

    with open_in(args.i1) as f1, open_in(args.i2) as f2, \
            open_out_gz(args.o1) as o1, open_out_gz(args.o2) as o2:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline().rstrip("\n")
            f1.readline()                       # '+'
            q1 = f1.readline().rstrip("\n")
            h2 = f2.readline()
            s2 = f2.readline().rstrip("\n")
            f2.readline()                       # '+'
            q2 = f2.readline().rstrip("\n")
            if not h2:
                break
            n_total += 1

            if args.anchor_linker:
                pos = find_linker(s1, linker, args.bc2_end, args.linker_max_shift)
                if pos < 0:
                    drop["no_linker"] += 1
                    continue
                shift = pos - args.bc2_end
                bc2 = s1[pos - bc2_len:pos]
                bc1 = s1[pos + Llen:pos + Llen + bc1_len]
                seq_start = args.seq_start + shift
            else:
                if len(s1) < args.bc1_end:
                    drop["r1_too_short_for_barcode"] += 1
                    continue
                bc2 = s1[args.bc2_start:args.bc2_end]
                bc1 = s1[args.bc1_start:args.bc1_end]
                seq_start = args.seq_start

            if len(bc2) != bc2_len or len(bc1) != bc1_len:
                drop["bad_barcode"] += 1
                continue

            barcode = (bc2 + bc1) if args.barcode_order == "BC2BC1" else (bc1 + bc2)

            g_seq = s1[seq_start:]
            g_qual = q1[seq_start:]
            if len(g_seq) < args.min_r1_len:
                drop["short_genomic_R1"] += 1
                continue

            base1, mate1 = name_base_and_mate(h1)
            _, mate2 = name_base_and_mate(h2)
            newname = "%s_CB=%s" % (base1, barcode)           # shared QNAME for the pair
            o1.write("@%s%s\n%s\n+\n%s\n" % (newname, mate1 or "/1", g_seq, g_qual))
            o2.write("@%s%s\n%s\n+\n%s\n" % (newname, mate2 or "/2", s2, q2))

            bc_counts[barcode] += 1
            n_kept += 1
            if args.preview_out and len(preview) < args.preview_n:
                preview.append((base1, bc2, bc1, barcode, g_seq[:20]))

    if args.barcode_counts_out:
        with open(args.barcode_counts_out, "w") as f:
            f.write("barcode\tcount\n")
            for bc, c in bc_counts.most_common():
                f.write("%s\t%d\n" % (bc, c))
    if args.barcode_whitelist_out:
        with open(args.barcode_whitelist_out, "w") as f:
            for bc, _ in bc_counts.most_common():
                f.write(bc + "\n")
    if args.drop_stats_out:
        with open(args.drop_stats_out, "w") as f:
            f.write("metric\tvalue\n")
            f.write("total_pairs\t%d\n" % n_total)
            f.write("kept_pairs\t%d\n" % n_kept)
            f.write("dropped_pairs\t%d\n" % (n_total - n_kept))
            for k, v in drop.most_common():
                f.write("dropped_%s\t%d\n" % (k, v))
            f.write("unique_barcodes\t%d\n" % len(bc_counts))
            f.write("kept_fraction\t%.4f\n" % ((n_kept / n_total) if n_total else 0.0))
    if args.preview_out:
        with open(args.preview_out, "w") as f:
            f.write("read\tBC2\tBC1\tbarcode\tgenomic_head20\n")
            for row in preview:
                f.write("\t".join(row) + "\n")

    sys.stderr.write("[make_trimmed_r1] total=%d kept=%d dropped=%d uniq_bc=%d\n"
                     % (n_total, n_kept, n_total - n_kept, len(bc_counts)))


if __name__ == "__main__":
    main()
