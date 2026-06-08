#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig1f_aggregate_cna.py

Replicate the slide-DNA-seq paper (Zhao et al., Nature 2022) Fig 1f/1j style:
aggregate ALL fragments into 1-Mb bins -> "Normalized copy number" along the genome,
baseline = 2 (diploid). Produces RAW vs GC+mappability-corrected panels so the user
can compare to the paper and see where the difference comes from.

Paper method (Suppl 3.2): 1Mb bins; select GC>0.35 & mappability>0.7; sequential
LOESS GC+mappability normalization; DIVIDE BY PAIRED BULK WGS (tissue-specific bias)
[WE LACK THIS]; modal-normalize (divide by mode of autosomal bins) x2 -> CN scaled to 2.
"""
import argparse, gzip, os, subprocess
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HG38 = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,
        "chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,
        "chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,
        "chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,
        "chr21":46709983,"chr22":50818468,"chrX":156040895}
CHROMS = ["chr%d"%i for i in range(1,23)]+["chrX"]
BIN = 1_000_000
PDAC = {"KRAS":("chr12",25.21),"CDKN2A":("chr9",21.97),"TP53":("chr17",7.67),
        "SMAD4":("chr18",51.03),"MYC":("chr8",127.74),"GATA6":("chr18",22.17)}

def build_bins():
    rows, off, gi = [], {}, 0
    for c in CHROMS:
        n=(HG38[c]+BIN-1)//BIN; off[c]=gi
        for b in range(n): rows.append((c,b*BIN,min((b+1)*BIN,HG38[c])))
        gi+=n
    return pd.DataFrame(rows,columns=["chr","start","end"]), off, gi

def bed_frac(path, off, n, skip_track):
    cov=np.zeros(n); op=gzip.open(path,"rt") if path.endswith(".gz") else open(path)
    with op as f:
        for line in f:
            if skip_track and line.startswith("track"): continue
            p=line.split("\t")
            if len(p)<3 or p[0] not in CHROMS: continue
            c,s,e=p[0],int(p[1]),int(p[2])
            for b in range(s//BIN,(e-1)//BIN+1):
                cov[off[c]+b]+=min(e,(b+1)*BIN)-max(s,b*BIN)
    return cov

def aggregate(frag, off, n):
    proc=subprocess.Popen(["gzip","-dc",frag],stdout=subprocess.PIPE)
    fr=pd.read_csv(proc.stdout,sep="\t",header=None,names=["chr","s","e","cb","q"],on_bad_lines="skip")
    proc.stdout.close(); proc.wait()
    fr=fr[fr.chr.isin(CHROMS)]
    mid=(fr.s.values+fr.e.values)//2
    coff=fr.chr.map(off).values
    cmax=fr.chr.map(lambda c:(HG38[c]+BIN-1)//BIN-1).values
    idx=(coff+np.minimum(mid//BIN,cmax)).astype(np.int64)
    return np.bincount(idx,minlength=n).astype(float)

def qcorrect(v, covar, keep, nq=30):
    out=v.copy().astype(float); vv=v[keep]; cc=covar[keep]
    q=pd.qcut(pd.Series(cc).rank(method="first"),nq,labels=False)
    med=pd.Series(vv).groupby(q.values).transform("median").values
    glob=np.median(vv[vv>0]) if (vv>0).any() else 1.0
    out[keep]=np.where(med>0, vv/med*glob, 0.0)
    return out

def modal_cn(v, keep):
    """divide by mode of autosomal bins (paper: 50 intervals), x2 -> CN scaled to 2."""
    x=v[keep]; x=x[x>0]
    hist,edges=np.histogram(x,bins=50)
    mode=(edges[np.argmax(hist)]+edges[np.argmax(hist)+1])/2
    cn=np.full(len(v),np.nan); cn[keep]=v[keep]/mode*2
    return cn

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--frag",required=True)
    ap.add_argument("--ref_dir",required=True)
    ap.add_argument("--sample",default="520_520")
    ap.add_argument("--control",default="",
                    help="optional matched bulk-WGS (or control) 1Mb coverage TSV "
                         "(cols chr,start,end,count) to DIVIDE BY, exactly like the paper "
                         "(visualize_coverage: norm_profiles = profiles./control). "
                         "Without it we replicate the paper's Fig S6b case (control=ones).")
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".",exist_ok=True)

    bins,off,n=build_bins()
    gcf=pd.read_csv(os.path.join(a.ref_dir,"bin_gc_1mb.tsv"),sep="\t")
    gc=(bins.chr+":"+(bins.start//BIN).astype(str)).map(
        dict(zip(gcf.chr+":"+(gcf.start//BIN).astype(str),gcf.gc))).fillna(0).values
    mapp=bed_frac(os.path.join(a.ref_dir,"k100.umap.bed.gz"),off,n,True)/BIN
    blk =bed_frac(os.path.join(a.ref_dir,"hg38-blacklist.v2.bed.gz"),off,n,False)/BIN

    agg=aggregate(a.frag,off,n)
    # paper-style bin selection: GC>0.35 & mappability>0.7 (+ our blacklist & coverage)
    cov_ok=agg>np.median(agg[agg>0])*0.25
    keep=(gc>0.35)&(mapp>0.7)&(blk<0.3)&cov_ok
    print(f"[{a.sample}] total fragments={int(agg.sum()):,}; bins kept={int(keep.sum())}/{n}")

    # RAW (no GC/map correction), modal-normalized to CN=2  -> shows GC/map waviness
    cn_raw=modal_cn(agg,keep)
    # GC then mappability quantile-correction (LOESS analogue), modal-normalized to 2
    c1=qcorrect(agg,gc,keep); c2=qcorrect(c1,mapp,keep)
    # OPTIONAL: divide by matched bulk-WGS control (paper's tissue-bias correction step)
    if a.control:
        ctrl=pd.read_csv(a.control,sep="\t",header=None,names=["chr","s","e","c"])
        cmap_={f"{r.chr}:{r.s//BIN}":r.c for r in ctrl.itertuples()}
        cvec=np.array([cmap_.get(f"{bins.chr[i]}:{bins.start[i]//BIN}",np.nan) for i in range(n)],float)
        cvec=qcorrect(cvec,gc,keep & ~np.isnan(cvec))     # bias-normalize control the same way
        cvec=qcorrect(cvec,mapp,keep & ~np.isnan(cvec))
        with np.errstate(divide="ignore",invalid="ignore"):
            c2=np.where((cvec>0)&keep, c2/cvec, np.nan)
        print(f"[{a.sample}] divided by control bulk-WGS: {a.control}")
    cn_cor=modal_cn(c2,keep)

    corr_raw=np.corrcoef(cn_raw[keep],gc[keep])[0,1]
    corr_cor=np.corrcoef(cn_cor[keep],gc[keep])[0,1]
    print(f"[{a.sample}] corr(CN,GC): raw={corr_raw:+.2f} -> corrected={corr_cor:+.2f}")

    # save table
    out_t=bins.copy(); out_t["keep"]=keep; out_t["cn_raw"]=cn_raw; out_t["cn_corrected"]=cn_cor
    out_t.to_csv(a.out+"_table.tsv",sep="\t",index=False)

    # chromosome layout (kept bins only, like the paper's concatenated genome axis)
    kb=bins[keep].reset_index(drop=True)
    cn_r=cn_raw[keep]; cn_c=cn_cor[keep]
    bnds,centers,order,pos=[],[],[],0
    for c in CHROMS:
        m=int((kb.chr.values==c).sum())
        if m: bnds.append(pos); centers.append(pos+m/2); order.append(c.replace("chr","")); pos+=m
    bnds.append(pos)
    gx={}
    for g,(gc_,mb) in PDAC.items():
        sel=kb.index[(kb.chr==gc_)&(kb.start<=mb*1e6)&(kb.end>mb*1e6)]
        if len(sel): gx[g]=sel[0]

    fig,axes=plt.subplots(2,1,figsize=(15,6),sharex=True)
    for ax,cn,ttl,cr in [(axes[0],cn_r,"RAW aggregate (no GC/mappability correction)",corr_raw),
                          (axes[1],cn_c,"GC + mappability corrected (our best, NO bulk-WGS division)",corr_cor)]:
        ax.scatter(range(len(cn)),cn,s=3,c=np.where(cn>2,"#c0392b","#2471a3"),alpha=0.6)
        ax.axhline(2,color="grey",lw=0.7); ax.set_ylim(0,4); ax.set_ylabel("Normalized\ncopy number")
        for b in bnds: ax.axvline(b-0.5,color="k",lw=0.3,alpha=0.4)
        ax.set_title(f"{a.sample}  —  {ttl}   [corr(CN,GC)={cr:+.2f}]",fontsize=10)
    for g,x in gx.items():
        axes[1].axvline(x,color="darkorange",lw=0.7,ls="--",alpha=0.8)
        axes[1].text(x,3.7,g,fontsize=7,rotation=90,va="top",color="darkorange")
    axes[1].set_xticks(centers); axes[1].set_xticklabels(order,fontsize=8)
    axes[1].set_xlabel("Genomic position (1 Mb bins, chr1-22,X)")
    fig.suptitle(f"{a.sample}: aggregate 1-Mb copy number (paper Fig 1f/1j style)",fontsize=12)
    plt.tight_layout(); plt.savefig(a.out+".png",dpi=200); plt.savefig(a.out+".pdf"); plt.close()
    print(f"[{a.sample}] wrote {a.out}.png / .pdf / _table.tsv")

if __name__=="__main__":
    main()
