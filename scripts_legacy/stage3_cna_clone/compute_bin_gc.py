#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compute_bin_gc.py <hg38.fa.gz> <out.tsv> [bin_size]
Compute GC fraction per fixed-size bin for chr1-22 + chrX (hg38).
Output columns: chr  start  end  gc  n_acgt
"""
import gzip
import sys

BIN = int(sys.argv[3]) if len(sys.argv) > 3 else 1_000_000
CHROMS = set("chr%d" % i for i in range(1, 23)) | {"chrX"}
fa, out = sys.argv[1], sys.argv[2]
rows = []


def flush(chrom, seq):
    if chrom not in CHROMS:
        return
    n = len(seq)
    for b in range(0, n, BIN):
        sub = seq[b:b + BIN]
        g = sub.count("G") + sub.count("C") + sub.count("g") + sub.count("c")
        a = sub.count("A") + sub.count("T") + sub.count("a") + sub.count("t")
        tot = g + a
        rows.append((chrom, b, min(b + BIN, n), (g / tot) if tot else 0.0, tot))


cur, buf = None, []
with gzip.open(fa, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if cur is not None:
                flush(cur, "".join(buf))
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        flush(cur, "".join(buf))

with open(out, "w") as o:
    o.write("chr\tstart\tend\tgc\tn_acgt\n")
    for r in rows:
        o.write("%s\t%d\t%d\t%.5f\t%d\n" % r)
print("[compute_bin_gc] wrote %d bins -> %s" % (len(rows), out))
