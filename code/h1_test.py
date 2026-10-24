
"""h1_test.py — preregistered H1 decision rule (notes/preregistration.md, 2026-10-23T08:34:58+01:00).
Usage: python code/h1_test.py <TEST_TARGET>   e.g. DP (discovery) or DA (held-out)
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy import stats
W="results/divergence"; T=sys.argv[1]; rng=np.random.default_rng(20261023)
df=pd.read_csv(f"{W}/gene_distances.csv")
need=["FH1","FH2","MH1","MH2",T]
ok=np.ones(len(df),bool)
for t in need:
    ok&=(df[f"{t}_cov"]>=0.8)&(df[f"{t}_ratio2"]<0.8)
d=df[ok&df.mh2_chr.notna()].copy()
def col(a,b): return f"d_{min(a,b)}_{max(a,b)}"
d["b"]=np.fmin(d["d_PAT_FH1"],d["d_PAT_FH2"])
d["t"]=np.fmin(d[col(T,"FH1")],d[col(T,"FH2")])
d["W"]=d[col("FH1","FH2")]
d["mat_floor"]=np.fmin(d["d_MAT_FH1"],d["d_MAT_FH2"])
d=d[d.b.notna()&d.t.notna()&d.W.notna()&d.mh2_start.notna()].sort_values(["mh2_chr","mh2_start"]).reset_index(drop=True)
n=len(d); out={"target":T,"n_genes":n}
q99=np.quantile(d.b,0.99); q95=np.quantile(d.b,0.95)
R=np.median(d.t)/np.median(d.b); f99=float((d.t>q99).mean()); f95=float((d.t>q95).mean())
out.update(dict(median_b=float(np.median(d.b)),median_t=float(np.median(d.t)),median_W=float(np.median(d.W)),R=float(R),q99_b=float(q99),q95_b=float(q95),f99=f99,f95=f95,
                median_mat_floor=float(np.nanmedian(d.mat_floor)),floor_ok=bool(np.nanmedian(d.mat_floor)<0.1*np.median(d.W))))
# block bootstrap (2-Mb blocks of genes along MH2)
d["blk"]=(d.mh2_chr.astype(int)*1000+(d.mh2_start//2_000_000)).astype(int)
blocks=d.blk.unique(); bi={b:np.where(d.blk.values==b)[0] for b in blocks}
boot=[]
for _ in range(1000):
    idx=np.concatenate([bi[b] for b in rng.choice(blocks,len(blocks))])
    tt=d.t.values[idx]; bb=d.b.values[idx]; qq=np.quantile(bb,0.99); q5=np.quantile(bb,0.95)
    boot.append((np.median(tt)/np.median(bb),(tt>qq).mean(),(tt>q5).mean()))
boot=np.array(boot)
out["R_CI"]=[float(x) for x in np.percentile(boot[:,0],[2.5,97.5])]; out["f99_CI"]=[float(x) for x in np.percentile(boot[:,1],[2.5,97.5])]; out["f95_CI"]=[float(x) for x in np.percentile(boot[:,2],[2.5,97.5])]
# spatial clustering of high-t genes: Moran's I at lag 1 along gene order within chromosomes vs circular shifts
hi=(d.t>q99).astype(float).values; chrs=d.mh2_chr.values
def moran(x):
    num=0.0; den=0.0; m=x.mean(); 
    for c in np.unique(chrs):
        xi=x[chrs==c]-m
        if len(xi)<3: continue
        num+=np.sum(xi[:-1]*xi[1:]); den+=np.sum(xi**2)
    return num/den if den>0 else np.nan
I_obs=moran(hi); null=[]
for _ in range(1000):
    xs=hi.copy()
    for c in np.unique(chrs):
        m_=chrs==c; xs[m_]=np.roll(hi[m_],rng.integers(1,max(2,m_.sum())))
    null.append(moran(xs))
null=np.array(null); out["moran_I"]=float(I_obs); out["moran_null99"]=float(np.nanpercentile(null,99)); out["moran_p"]=float((np.sum(null>=I_obs)+1)/(len(null)+1))
# secondary: d(T, DA) bimodality if available
other="DA" if T=="DP" else "DP"
c2=col(T,other)
if c2 in df.columns:
    x=df.loc[ok, c2].dropna().values
    try:
        import diptest; dip,p=diptest.diptest(x); out["dip_TvsOther"]=float(dip); out["dip_p"]=float(p)
    except Exception as e: out["dip_err"]=str(e)
    out["n_TvsOther"]=int(len(x)); out["median_TvsOther"]=float(np.median(x)) if len(x) else None
    out["frac_TvsOther_gt_q99b"]=float((x>q99).mean()) if len(x) else None
# scale control
if col(T,"UU") in df.columns or "d_FH1_UU" in df.columns:
    uu=df.loc[ok,"d_FH1_UU"].dropna(); out["median_d_FH1_UU"]=float(uu.median()); out["scale_ok"]=bool(uu.median()>=3*np.median(d.b))
# verdict per preregistered rule
same=(0.67<=R<=1.5) and f99<=0.05 and not (out.get("dip_p",1)<0.01 and out.get("frac_TvsOther_gt_q99b",0)>0.05)
allo1=f99>=0.20 and I_obs>out["moran_null99"] and out.get("dip_p",1)<0.01
allo2=R>=2 and f99>=0.50
out["verdict"]="same-lineage tetraploid" if same else ("allotetraploid (one subgenome shared)" if allo1 else ("allotetraploid (both subgenomes divergent)" if allo2 else "inconclusive"))
# per-chromosome R
out["per_chr_R"]={int(c):float(np.median(g.t)/np.median(g.b)) for c,g in d.groupby("mh2_chr")}
out["per_chr_f99"]={int(c):float((g.t>q99).mean()) for c,g in d.groupby("mh2_chr")}
d.to_csv(f"{W}/h1_genes_{T}.csv",index=False)
json.dump(out,open(f"{W}/h1_result_{T}.json","w"),indent=1)
print(json.dumps({k:v for k,v in out.items() if not k.startswith("per_chr")},indent=1))
