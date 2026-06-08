#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage3_per_spot_cnv.py

Per-spot copy-number (CNV) for spatial DNA-seq (DBiT 50x50, hg38).

A single 50um spot is too sparse (~1 fragment / 1Mb bin) for direct per-bin copy
number, so we use SPATIAL kNN SMOOTHING (paper's approach): each spot's CNV is
computed from itself + its k nearest spatial neighbours. The effective resolution
= the smoothing neighbourhood (tune via --spatial_win / --bin_size).

Pipeline: spot x bin counts -> bin QC -> per-spot GC + mappability correction
-> spatial smoothing (sum over neighbours) -> library norm -> per-bin pseudo-normal
reference -> per-spot modal anchor x2 -> copy number. Optional per-chrom median
segmentation.

Outputs:
  per_spot_cnv_matrix.tsv.gz   rows = spots (x_id,y_id,total_frags), cols = 1Mb bins -> copy number
  bin_info.tsv                 bin -> chr,start,end (column order of the matrix)
  per_spot_cnv_burden.tsv      per spot: fraction of genome with |CN-2|>cna_thresh
  spatial_cn_<chr>.png         copy number painted on the 50x50 grid, per requested chromosome
  example_spot_profiles.png    genome-wide CN for a few example spots
"""
import argparse, gzip, os, subprocess
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.ndimage import uniform_filter, median_filter

HG38 = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,
        "chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,
        "chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,
        "chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,
        "chr21":46709983,"chr22":50818468,"chrX":156040895}
CHROMS = ["chr%d"%i for i in range(1,23)]+["chrX"]
AUTO = ["chr%d"%i for i in range(1,23)]
GRID, BIN = 50, 1_000_000


def build_bins(bs):
    rows, off, gi = [], {}, 0
    for c in CHROMS:
        n=(HG38[c]+bs-1)//bs; off[c]=gi
        for b in range(n): rows.append((c,b*bs,min((b+1)*bs,HG38[c])))
        gi+=n
    return pd.DataFrame(rows,columns=["chr","start","end"]), off, gi


def bed_frac(path, off, n, bs, skip):
    cov=np.zeros(n); op=gzip.open(path,"rt") if path.endswith(".gz") else open(path)
    with op as f:
        for line in f:
            if skip and line.startswith("track"): continue
            p=line.split("\t")
            if len(p)<3 or p[0] not in CHROMS: continue
            c,s,e=p[0],int(p[1]),int(p[2])
            for b in range(s//bs,(e-1)//bs+1):
                cov[off[c]+b]+=min(e,(b+1)*bs)-max(s,b*bs)
    return cov


def load_counts(frag, matched, off, n, bs):
    m=pd.read_csv(matched,sep="\t",usecols=["barcode_obs","x_id","y_id","full_status"])
    m=m[m.full_status=="mapped"].dropna(subset=["x_id","y_id"])
    ox=dict(zip(m.barcode_obs,m.x_id.astype(int))); oy=dict(zip(m.barcode_obs,m.y_id.astype(int)))
    pr=subprocess.Popen(["gzip","-dc",frag],stdout=subprocess.PIPE)
    fr=pd.read_csv(pr.stdout,sep="\t",header=None,names=["chr","s","e","cb","q"],on_bad_lines="skip")
    pr.stdout.close(); pr.wait()
    fr=fr[fr.chr.isin(CHROMS)]; fr["x"]=fr.cb.map(ox); fr["y"]=fr.cb.map(oy)
    fr=fr.dropna(subset=["x","y"]); fr["x"]=fr.x.astype(int); fr["y"]=fr.y.astype(int)
    fr=fr[(fr.x>=1)&(fr.x<=GRID)&(fr.y>=1)&(fr.y<=GRID)]
    mid=(fr.s.values+fr.e.values)//2; coff=fr.chr.map(off).values
    cmax=fr.chr.map(lambda c:(HG38[c]+bs-1)//bs-1).values
    bidx=(coff+np.minimum(mid//bs,cmax)).astype(np.int64)
    sidx=((fr.x.values-1)*GRID+(fr.y.values-1)).astype(np.int64)
    return np.bincount(sidx*n+bidx,minlength=GRID*GRID*n).reshape(GRID*GRID,n)


def qcorrect_rows(C, covar, n_q=20):
    """per-spot flatten of counts vs a bin-level covariate (GC or mappability)."""
    ranks=np.argsort(np.argsort(covar)); grp=np.minimum(ranks*n_q//len(covar),n_q-1)
    oh=np.zeros((len(covar),n_q)); oh[np.arange(len(covar)),grp]=1.0
    gc=oh.sum(0); gs=C@oh
    exp=np.divide(gs,gc[None,:],out=np.zeros_like(gs),where=gc[None,:]>0)[:,grp]
    rm=C.mean(1,keepdims=True)
    return np.divide(C,exp,out=np.zeros_like(C),where=exp>0)*rm


def spatial_sum(mat, present, win):
    nb=mat.shape[1]; data=np.zeros((GRID,GRID,nb)); mask=np.zeros((GRID,GRID))
    for s in range(GRID*GRID):
        if present[s]: data[s//GRID,s%GRID]=mat[s]; mask[s//GRID,s%GRID]=1
    num=uniform_filter(data*mask[:,:,None],size=(win,win,1),mode="constant")
    return num.reshape(GRID*GRID,nb)   # neighbourhood SUM of corrected counts


def resolve_sample(stage2_dir, s):
    """find the Stage-2 dir for sample s (auto-detect <s>/, <s>_xy_map/, <s>_out/)."""
    for cand in (s, s+"_xy_map", s+"_out"):
        d=os.path.join(stage2_dir,cand)
        if os.path.isfile(os.path.join(d,"top50x50_fragments_with_cb.tsv.gz")): return d
    return os.path.join(stage2_dir,s)


def load_normal_mask(path):
    """read a normal/reference spot list (cols x_id,y_id) -> boolean mask over GRID*GRID, or None."""
    if not path or not os.path.isfile(path): return None
    df=pd.read_csv(path,sep="\t")
    if "x_id" not in df or "y_id" not in df: return None
    sidx=((df.x_id.astype(int)-1)*GRID+(df.y_id.astype(int)-1)).values
    sidx=sidx[(sidx>=0)&(sidx<GRID*GRID)]
    m=np.zeros(GRID*GRID,bool); m[sidx]=True
    return m


def run_sample(frag, matched, sample, outdir, ref_dir, bs,
               spatial_win, genome_win, cna_thresh, cn_chroms,
               gc, mapp, blk, bins, off, n, normal_mask=None):
    os.makedirs(outdir,exist_ok=True)
    counts=load_counts(frag,matched,off,n,bs)
    tot=counts.sum(1).astype(float); present=tot>0
    bin_mean=counts[present].mean(0)
    keep=(gc>0.35)&(mapp>0.5)&(blk<0.3)&(bin_mean>0.2*np.median(bin_mean[bin_mean>0]))
    kb=bins[keep].reset_index(drop=True); kb.to_csv(os.path.join(outdir,"bin_info.tsv"),sep="\t",index=False)
    print("[%s] spots=%d bins kept=%d/%d  (spatial_win=%d ~ kNN%d)"%(sample,int(present.sum()),int(keep.sum()),n,spatial_win,spatial_win**2))

    C=counts[:,keep].astype(float); G=gc[keep]; Mp=mapp[keep]
    # per-spot GC + mappability correction
    C=qcorrect_rows(C,G); C=qcorrect_rows(C,Mp)
    # SPATIAL SMOOTHING: neighbourhood sum (this is what makes per-spot CNV usable)
    Cs=spatial_sum(C,present,spatial_win)
    # library normalize each (smoothed) spot
    t2=Cs.sum(1); med=np.median(t2[present]); libn=np.zeros_like(Cs)
    libn[present]=Cs[present]/t2[present,None]*med
    # per-bin pseudo-normal reference -> relative.
    # if normal_mask given (histology stroma / data-driven flattest spots), build the
    # reference ONLY from those spots; else fall back to cross-spot median (all spots).
    ref_spots=present
    if normal_mask is not None:
        nm=present & normal_mask
        if nm.sum()>=20: ref_spots=nm; print("[%s] normal-anchored reference: %d spots"%(sample,int(nm.sum())))
        else: print("[%s] WARN normal_spots too few (%d<20); using all-spot median"%(sample,int(nm.sum())))
    ref=np.median(libn[ref_spots],0); ref[ref<=0]=np.nan
    rel=np.full_like(libn,np.nan); rel[present]=libn[present]/ref
    # per-spot modal anchor (autosomal median -> 1) x2 -> copy number
    autos=kb.chr.isin(AUTO).values
    base=np.nanmedian(rel[:,autos],axis=1,keepdims=True); base[base<=0]=np.nan
    CN=rel/base*2.0
    # genomic smoothing per chrom
    for c in CHROMS:
        idx=np.where(kb.chr.values==c)[0]
        if len(idx)>=genome_win:
            CN[np.ix_(present,idx)]=median_filter(CN[np.ix_(present,idx)],size=(1,genome_win),mode="nearest")
    CN=np.clip(CN,0,6)

    # ---- outputs ----
    xs=np.repeat(np.arange(1,GRID+1),GRID); ys=np.tile(np.arange(1,GRID+1),GRID)
    cn_df=pd.DataFrame(CN,columns=[f"{r.chr}:{r.start}" for r in kb.itertuples()])
    cn_df.insert(0,"total_frags",tot.astype(int)); cn_df.insert(0,"y_id",ys); cn_df.insert(0,"x_id",xs)
    cn_df[present].to_csv(os.path.join(outdir,"per_spot_cnv_matrix.tsv.gz"),sep="\t",index=False,compression="gzip")
    burden=np.full(GRID*GRID,np.nan)
    burden[present]=(np.abs(CN[present]-2)>cna_thresh).mean(1)
    pd.DataFrame({"x_id":xs,"y_id":ys,"total_frags":tot.astype(int),"cnv_burden":burden}).to_csv(
        os.path.join(outdir,"per_spot_cnv_burden.tsv"),sep="\t",index=False)

    # spatial CN maps per requested chromosome
    for c in [x.strip() for x in cn_chroms.split(",") if x.strip()]:
        idx=np.where(kb.chr.values==c)[0]
        if len(idx)==0: continue
        cn_c=np.nanmean(CN[:,idx],axis=1)
        grid=np.full((GRID,GRID),np.nan)
        for s in range(GRID*GRID):
            if present[s]: grid[s//GRID,s%GRID]=cn_c[s]
        plt.figure(figsize=(7,6))
        im=plt.imshow(grid,origin="lower",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=1,vcenter=2,vmax=3),interpolation="nearest")
        plt.colorbar(im,label="copy number"); plt.xlabel("Y"); plt.ylabel("X")
        plt.title("%s: per-spot copy number — %s"%(sample,c))
        plt.tight_layout(); plt.savefig(os.path.join(outdir,"spatial_cn_%s.png"%c),dpi=200); plt.close()

    # example spot genome-wide profiles
    pres_idx=np.where(present)[0]
    sel=pres_idx[np.linspace(0,len(pres_idx)-1,4).astype(int)]
    bnds,centers,order,pos=[],[],[],0
    for c in CHROMS:
        m=int((kb.chr.values==c).sum())
        if m: bnds.append(pos); centers.append(pos+m/2); order.append(c.replace("chr","")); pos+=m
    fig,axes=plt.subplots(len(sel),1,figsize=(14,2*len(sel)),squeeze=False)
    for ax,s in zip(axes[:,0],sel):
        ax.plot(CN[s],lw=0.5,color="#333"); ax.axhline(2,color="grey",lw=0.5); ax.set_ylim(0,5)
        for b in bnds: ax.axvline(b-0.5,color="k",lw=0.3,alpha=0.4)
        ax.set_xticks(centers); ax.set_xticklabels(order,fontsize=7)
        ax.set_ylabel("spot(%d,%d)\nCN"%(s//GRID+1,s%GRID+1),fontsize=8)
    axes[0,0].set_title("%s: example per-spot genome-wide copy number (spatially smoothed)"%sample)
    plt.tight_layout(); plt.savefig(os.path.join(outdir,"example_spot_profiles.png"),dpi=200); plt.close()
    print("[%s] DONE -> %s  (per_spot_cnv_matrix.tsv.gz: %d spots x %d bins)"%(sample,outdir,int(present.sum()),int(keep.sum())))


def main():
    ap=argparse.ArgumentParser(description="Per-spot CNV for spatial DNA-seq (DBiT 50x50, hg38).")
    ap.add_argument("--stage2_dir",required=True,help="dir holding Stage-2 outputs (auto-detects <sample>/, <sample>_xy_map/, <sample>_out/)")
    ap.add_argument("--samples",nargs="+",default=["520_520"])
    ap.add_argument("--ref_dir",required=True)
    ap.add_argument("--outroot",required=True)
    ap.add_argument("--bin_size",type=int,default=1_000_000)
    ap.add_argument("--spatial_win",type=int,default=7,help="NxN neighbourhood (7=~kNN50); larger=cleaner CNV, coarser spatial")
    ap.add_argument("--genome_win",type=int,default=3)
    ap.add_argument("--cna_thresh",type=float,default=0.5)
    ap.add_argument("--cn_chroms",default="chr8,chr18",help="chromosomes to paint as spatial CN maps")
    ap.add_argument("--normal_spots",default=None,
                    help="TSV (cols x_id,y_id) of reference/normal spots for the pseudo-normal; "
                         "use {sample} placeholder for per-sample files. If omitted, uses all-spot median.")
    a=ap.parse_args(); bs=a.bin_size

    bins,off,n=build_bins(bs)
    gcf=pd.read_csv(os.path.join(a.ref_dir,"bin_gc_1mb.tsv"),sep="\t")
    gc=(bins.chr+":"+(bins.start//BIN).astype(str)).map(
        dict(zip(gcf.chr+":"+(gcf.start//BIN).astype(str),gcf.gc))).fillna(0).values
    mapp=bed_frac(os.path.join(a.ref_dir,"k100.umap.bed.gz"),off,n,bs,True)/bs
    blk =bed_frac(os.path.join(a.ref_dir,"hg38-blacklist.v2.bed.gz"),off,n,bs,False)/bs

    for s in a.samples:
        sd=resolve_sample(a.stage2_dir,s)
        frag=os.path.join(sd,"top50x50_fragments_with_cb.tsv.gz")
        matched=os.path.join(sd,"matched_spot_barcodes.tsv")
        if not os.path.isfile(frag):
            print("[%s] SKIP — not found under %s"%(s,a.stage2_dir)); continue
        npath=a.normal_spots.format(sample=s) if (a.normal_spots and "{sample}" in a.normal_spots) else a.normal_spots
        nmask=load_normal_mask(npath)
        try:
            run_sample(frag,matched,s,os.path.join(a.outroot,s),a.ref_dir,bs,
                       a.spatial_win,a.genome_win,a.cna_thresh,a.cn_chroms,
                       gc,mapp,blk,bins,off,n,normal_mask=nmask)
        except Exception as ex:
            print("[%s] FAILED: %s"%(s,ex))


if __name__=="__main__":
    main()
