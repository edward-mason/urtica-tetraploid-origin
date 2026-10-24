
"""reviewerA_tier1.py — independent recomputation of the H1 statistics (Reviewer A).
Written without reference to h1_test.py's internals beyond the stated filters.
"""
import sys, json, numpy as np, pandas as pd
GD=sys.argv[1]; OUT=sys.argv[2]
df=pd.read_csv(GD)
# --- independently re-derive parental assignment from Hirabayashi 2026 Table 2 (maternal = MH1 on 1,2,4,7-13; MH2 on 3,5,6)
MAT={c:("MH2" if c in (3,5,6) else "MH1") for c in range(1,14)}
def dcol(a,b): a,b=sorted([a,b]); return f"d_{a}_{b}"
chrom=df.mh2_chr
pat=[("MH1" if MAT[int(c)]=="MH2" else "MH2") if pd.notna(c) else None for c in chrom]
my_pat_fh1=np.array([df.loc[i,dcol(p,"FH1")] if p else np.nan for i,p in enumerate(pat)])
my_pat_fh2=np.array([df.loc[i,dcol(p,"FH2")] if p else np.nan for i,p in enumerate(pat)])
chk1=np.nanmax(np.abs(my_pat_fh1-df.d_PAT_FH1.values)); chk2=np.nanmax(np.abs(my_pat_fh2-df.d_PAT_FH2.values))
res={"parental_assignment_check_maxabsdiff":[float(chk1),float(chk2)]}
def analyse(T, agg="min", mode="standard", min_cols=None, excl_chr=None, nboot=2000, seed=7):
    need=["FH1","FH2","MH1","MH2",T]
    ok=np.ones(len(df),bool)
    for t in need: ok&=(df[f"{t}_cov"].values>=0.8)&(df[f"{t}_ratio2"].values<0.8)
    ok&=df.mh2_chr.notna().values
    d=df[ok].copy()
    if mode=="standard":
        f=np.fmin if agg=="min" else np.fmax
        b=f(my_pat_fh1[ok],my_pat_fh2[ok]); t=f(d[dcol(T,"FH1")].values,d[dcol(T,"FH2")].values)
        ncol=np.fmin(np.fmin(d[f"n_FH1_FH2"].values,d["n_"+dcol(T,"FH1")[2:]].values),d["n_"+dcol(T,"FH2")[2:]].values)
    else:  # (b) t' = d(T, PAT), b' = d(PAT, FH1)
        t=np.array([d.loc[i,dcol(T,p)] for i,p in zip(d.index,[pat[i] for i in d.index])])
        b=my_pat_fh1[ok]
        ncol=d["n_"+dcol(T,"MH1")[2:]].values
    keep=~np.isnan(b)&~np.isnan(t)&d.mh2_start.notna().values
    if min_cols: keep&=(ncol>=min_cols)
    if excl_chr: keep&=(d.mh2_chr.values!=excl_chr)
    b=b[keep]; t=t[keep]; d=d[keep]
    q99=np.quantile(b,0.99); R=np.median(t)/np.median(b); f99=float((t>q99).mean())
    # bootstrap of my own design: resample chromosome-arm-scale blocks (5-Mb windows, MH2 coords) AND a plain gene bootstrap
    rng=np.random.default_rng(seed)
    blk=(d.mh2_chr.values.astype(int)*10000+(d.mh2_start.values//5_000_000)).astype(int)
    ub,inv=np.unique(blk,return_inverse=True); members=[np.where(inv==k)[0] for k in range(len(ub))]
    bR=[];bf=[]
    for _ in range(nboot):
        idx=np.concatenate([members[k] for k in rng.integers(0,len(ub),len(ub))])
        bb=b[idx];tt=t[idx]; bR.append(np.median(tt)/np.median(bb)); bf.append((tt>np.quantile(bb,0.99)).mean())
    gR=[];gf=[]
    for _ in range(nboot):
        idx=rng.integers(0,len(b),len(b)); bb=b[idx];tt=t[idx]; gR.append(np.median(tt)/np.median(bb)); gf.append((tt>np.quantile(bb,0.99)).mean())
    return dict(n=int(len(b)),median_b=float(np.median(b)),median_t=float(np.median(t)),q99_b=float(q99),R=float(R),f99=f99,
                R_CI_block5Mb=[float(x) for x in np.percentile(bR,[2.5,97.5])],f99_CI_block5Mb=[float(x) for x in np.percentile(bf,[2.5,97.5])],
                R_CI_gene=[float(x) for x in np.percentile(gR,[2.5,97.5])],f99_CI_gene=[float(x) for x in np.percentile(gf,[2.5,97.5])],
                n_blocks=int(len(ub)))
for T in ["DP","DA"]:
    res[f"{T}_main"]=analyse(T)
    res[f"{T}_a_max"]=analyse(T,agg="max")
    res[f"{T}_b_PATbased"]=analyse(T,mode="pat")
    res[f"{T}_c_min1000cols"]=analyse(T,min_cols=1000)
    res[f"{T}_d_noChr8"]=analyse(T,excl_chr=8)
json.dump(res,open(OUT,"w"),indent=1)
for k,v in res.items(): print(k, v if not isinstance(v,dict) else {kk:v[kk] for kk in ("n","R","f99","R_CI_block5Mb","f99_CI_block5Mb")})
