
# -*- coding: utf-8 -*-
"""h1_test_v2.py — H1 statistics after independent review (2026-10-24).
Changes vs h1_test.py (all declared in the manuscript): (1) spatial clustering by within-chromosome permutation of the
high-t indicator (adjacent-pair count), replacing the degenerate circular-shift Moran's I; (2) windowed paired t-b scan
with permutation null; (3) paralog filter: genes whose target-vs-U.urens distance exceeds the FH1-vs-U.urens distance by
> 5 percentage points are flagged and the statistics are reported with and without them; (4) chromosome-level bootstrap
CIs in addition to 2-Mb blocks; (5) dip test executed in-script; (6) (delta, phi) power model.
Usage: python code/h1_test_v2.py DP|DA
"""
import sys, os, json
import numpy as np, pandas as pd
import diptest
W="results/divergence"; T=sys.argv[1]; rng=np.random.default_rng(20261023)
df=pd.read_csv(f"{W}/gene_distances.csv")
ok=np.ones(len(df),bool)
for t in ["FH1","FH2","MH1","MH2",T]: ok&=(df[f"{t}_cov"]>=0.8)&(df[f"{t}_ratio2"]<0.8)
d=df[ok&df.mh2_chr.notna()&df.mh2_start.notna()].copy()
def col(a,b): return f"d_{min(a,b)}_{max(a,b)}"
d["b"]=np.fmin(d.d_PAT_FH1,d.d_PAT_FH2); d["t"]=np.fmin(d[col(T,"FH1")],d[col(T,"FH2")]); d["W"]=d.d_FH1_FH2
d["dTUU"]=d[col(T,"UU")]; d["dFUU"]=d["d_FH1_UU"]
d=d[d.b.notna()&d.t.notna()&d.W.notna()].sort_values(["mh2_chr","mh2_start"]).reset_index(drop=True)
d["paralog_flag"]=(d.dTUU-d.dFUU)>0.05   # target locus far more divergent from U. urens than the diploid ortholog is
out={"target":T,"n_genes":int(len(d)),"n_paralog_flag":int(d.paralog_flag.sum())}
def stats(dd):
    q99=np.quantile(dd.b,.99); return dict(n=int(len(dd)),median_b=float(np.median(dd.b)),median_t=float(np.median(dd.t)),R=float(np.median(dd.t)/np.median(dd.b)),q99_b=float(q99),f99=float((dd.t>q99).mean()),f95=float((dd.t>np.quantile(dd.b,.95)).mean()),
                                        median_paired_diff=float(np.median(dd.t-dd.b)),mean_paired_diff=float(np.mean(dd.t-dd.b)))
out["all"]=stats(d); dc=d[~d.paralog_flag].copy(); out["paralog_filtered"]=stats(dc)
# bootstrap CIs: 2-Mb blocks and whole chromosomes (on paralog-filtered set)
def boot(dd,key,B=1000):
    groups=dd[key].values; ug=np.unique(groups); idx={g:np.where(groups==g)[0] for g in ug}; res=[]
    for _ in range(B):
        ii=np.concatenate([idx[g] for g in rng.choice(ug,len(ug))]); tt=dd.t.values[ii]; bb=dd.b.values[ii]; q=np.quantile(bb,.99)
        res.append((np.median(tt)/np.median(bb),(tt>q).mean()))
    res=np.array(res); return dict(R_CI=[float(x) for x in np.percentile(res[:,0],[2.5,97.5])],f99_CI=[float(x) for x in np.percentile(res[:,1],[2.5,97.5])])
dc["blk"]=(dc.mh2_chr.astype(int)*1000+(dc.mh2_start//2_000_000)).astype(int)
out["CI_2Mb_blocks"]=boot(dc,"blk"); out["CI_chromosomes"]=boot(dc,"mh2_chr")
# clustering test: adjacent pairs of high-t genes vs within-chromosome permutation
q99=out["paralog_filtered"]["q99_b"]; hi=(dc.t>q99).values.astype(int); chrs=dc.mh2_chr.values
def adj_pairs(x):
    s=0
    for c in np.unique(chrs):
        xi=x[chrs==c]; s+=int(np.sum(xi[:-1]*xi[1:]))
    return s
obs=adj_pairs(hi); null=[]
for _ in range(2000):
    xs=hi.copy()
    for c in np.unique(chrs):
        m=chrs==c; xs[m]=rng.permutation(hi[m])
    null.append(adj_pairs(xs))
null=np.array(null); out["clustering"]=dict(n_high=int(hi.sum()),adjacent_pairs_obs=obs,adjacent_pairs_null_mean=float(null.mean()),p=float((np.sum(null>=obs)+1)/(len(null)+1)))
# windowed paired scan: rolling 50-gene median of (t-b), null by permuting genes within chromosome
diff=(dc.t-dc.b).values; win=50
def roll_max(x):
    mx=-1; where=None
    for c in np.unique(chrs):
        m=np.where(chrs==c)[0]
        if len(m)<win: continue
        xi=x[m]; cs=np.convolve(xi,np.ones(win)/win,mode="valid")
        # use rolling median for robustness
        med=pd.Series(xi).rolling(win).median().values[win-1:]
        j=int(np.nanargmax(med))
        if med[j]>mx: mx=float(med[j]); where=(int(c),int(dc.mh2_start.values[m[j]]),int(dc.mh2_start.values[m[j+win-1]]))
    return mx,where
obs_mx,obs_where=roll_max(diff); nulls=[]
for _ in range(300):
    xs=diff.copy()
    for c in np.unique(chrs):
        m=chrs==c; xs[m]=rng.permutation(diff[m])
    nulls.append(roll_max(xs)[0])
thr=float(np.percentile(nulls,99))
# list all windows exceeding the null 99th percentile (positive = t elevated relative to b)
blocks=[]
for c in np.unique(chrs):
    m=np.where(chrs==c)[0]
    if len(m)<win: continue
    med=pd.Series(diff[m]).rolling(win).median().values
    for j in range(win-1,len(m)):
        if med[j]>thr:
            g=dc.iloc[m[j-win+1:j+1]]
            blocks.append(dict(chr=int(c),start=int(g.mh2_start.min()),end=int(g.mh2_start.max()),med_t=float(np.median(g.t)),med_b=float(np.median(g.b)),med_W=float(np.median(g.W)),frac_b0=float((g.b==0).mean())))
# merge overlapping windows per chromosome
merged=[]
for bk in sorted(blocks,key=lambda x:(x["chr"],x["start"])):
    if merged and merged[-1]["chr"]==bk["chr"] and bk["start"]<=merged[-1]["end"]:
        merged[-1]["end"]=max(merged[-1]["end"],bk["end"]); merged[-1]["windows"]+=1
    else: bk=dict(bk); bk["windows"]=1; merged.append(bk)
out["paired_scan"]=dict(window=win,null99_rolling_median_diff=thr,max_obs=obs_mx,max_where=obs_where,blocks=merged,
                        max_windowed_median_t=float(max(pd.Series(dc.t.values[chrs==c]).rolling(win).median().max() for c in np.unique(chrs) if (chrs==c).sum()>=win)))
# dip tests (raw and log) for d(T,other), W, b
other="DA" if T=="DP" else "DP"; x=df.loc[ok,col(T,other)].dropna().values
dips={}
for name,arr in [("d_T_other",x),("W",dc.W.values),("b",dc.b.values)]:
    dp,pp=diptest.diptest(arr); xl=np.log10(arr[arr>0]); dl,pl=diptest.diptest(xl)
    dips[name]=dict(n=int(len(arr)),dip_raw=float(dp),p_raw=float(pp),dip_log=float(dl),p_log=float(pl),median=float(np.median(arr)),frac_gt_q99b=float((arr>q99).mean()))
out["dip"]=dips
# power model: fraction phi of genes receive extra divergence delta on top of a b-like baseline
sim={}; bvals=dc.b.values
for delta in [0.0025,0.005,0.01,0.015,0.02,0.03,0.05]:
    for phi in [0.1,0.25,0.5,0.75,1.0]:
        Rs=[];Fs=[]
        for _ in range(20):
            tt=rng.choice(bvals,len(bvals)); mask=rng.random(len(bvals))<phi; tt=tt+mask*delta
            Rs.append(np.median(tt)/np.median(bvals)); Fs.append((tt>q99).mean())
        sim[f"delta={delta},phi={phi}"]=dict(R=float(np.mean(Rs)),f99=float(np.mean(Fs)))
out["power_model"]=sim
# excluded region: (delta,phi) combos whose expected R or f99 exceeds the observed upper CI
Ru=out["CI_chromosomes"]["R_CI"][1]; Fu=out["CI_chromosomes"]["f99_CI"][1]
out["excluded_combos"]=[k for k,v in sim.items() if v["R"]>Ru or v["f99"]>Fu]
# symmetric-design subset (chromosomes 3,5,6 where PAT=MH1 is not the query)
sym=dc[dc.mh2_chr.isin([3,5,6])]; out["symmetric_chr356"]=stats(sym); out["other_chr"]=stats(dc[~dc.mh2_chr.isin([3,5,6])])
# verdict (preregistered rule applied to paralog-filtered set, dip on log scale with diploid baselines as reference)
pf=out["paralog_filtered"]
same=(0.67<=pf["R"]<=1.5) and pf["f99"]<=0.05 and not (dips["d_T_other"]["p_log"]<0.01 and dips["d_T_other"]["frac_gt_q99b"]>0.05)
allo1=pf["f99"]>=0.20 and out["clustering"]["p"]<0.01 and dips["d_T_other"]["p_log"]<0.01
allo2=pf["R"]>=2 and pf["f99"]>=0.5
out["verdict"]="same-lineage tetraploid" if same else ("allotetraploid (one subgenome shared)" if allo1 else ("allotetraploid (both subgenomes divergent)" if allo2 else "inconclusive"))
dc.to_csv(f"{W}/h1v2_genes_{T}.csv",index=False); json.dump(out,open(f"{W}/h1v2_result_{T}.json","w"),indent=1)
print(json.dumps({k:v for k,v in out.items() if k not in ("power_model",)},indent=1)[:6000])
