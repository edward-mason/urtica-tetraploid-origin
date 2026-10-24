
"""h3_test.py — preregistered H3: diversity landscape vs predicted centromere class.
Usage: python code/h3_test.py female   (discovery: pi from FH1,FH2,PAT; positions FH1; Table S8 H1)
       python code/h3_test.py male     (held-out: pi = d(MH1,MH2); positions MH2; male centromere predictions file)
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy import stats
from itertools import combinations
W="results/divergence"; mode=sys.argv[1]; rng=np.random.default_rng(20261023)
df=pd.read_csv(f"{W}/gene_distances.csv")
cen=pd.read_csv("data/annot/centromere_points.csv")
ACRO={8,9,11,12,13}; META={10}
ok=(df.FH1_cov>=0.8)&(df.FH1_ratio2<0.8)&(df.FH2_cov>=0.8)&(df.FH2_ratio2<0.8)&(df.MH1_cov>=0.8)&(df.MH1_ratio2<0.8)&(df.MH2_cov>=0.8)&(df.MH2_ratio2<0.8)
d=df[ok].copy()
if mode=="female":
    d["pi"]=d[["d_FH1_FH2","d_PAT_FH1","d_PAT_FH2"]].mean(axis=1); d["chr"]=d.FH1_chr; d["pos"]=d.FH1_start; cp=cen[cen.hap=="H1"].set_index("chr")
else:
    d["pi"]=d["d_MH1_MH2"]; d["chr"]=d.MH2_chr; d["pos"]=d.MH2_start
    cp=pd.read_csv("data/annot/centromere_points_male.csv").set_index("chr")
d=d[d.pi.notna()&d.chr.notna()&d.pos.notna()].copy(); d["chr"]=d.chr.astype(int)
# chromosome lengths from the positions (max gene position as proxy) — replaced by assembly lengths if available
lens=json.load(open("data/annot/chrom_lengths.json")) if os.path.exists("data/annot/chrom_lengths.json") else {}
key="FH1" if mode=="female" else "MH2"
res=[]
for c,g in d.groupby("chr"):
    if c not in cp.index: continue
    cpt=cp.loc[c,"cen_point"]; L=lens.get(key,{}).get(str(c), g.pos.max())
    g=g.sort_values("pos").copy()
    g["dist_cen"]=np.abs(g.pos-cpt)/1e6; g["dist_end"]=np.minimum(g.pos,L-g.pos)/1e6
    rho,p=stats.spearmanr(g.dist_cen,g.pi); rho_end,p_end=stats.spearmanr(g.dist_end,g.pi)
    # circular-shift null for rho (keeps spatial autocorrelation of pi)
    pi=g.pi.values; null=[]
    for _ in range(1000):
        s=np.roll(pi,rng.integers(1,len(pi))); null.append(stats.spearmanr(g.dist_cen,s)[0])
    null=np.array(null); p_shift=(np.sum(null>=rho)+1)/(len(null)+1)
    k=max(5,int(0.2*len(g))); near=g.nsmallest(k,"dist_cen"); far=g.drop(near.index)
    Tr=np.median(near.pi)/np.median(far.pi)
    # preregistered control: remove genes within 2 Mb of the predicted centromere point (gene-density circularity), recompute T
    g2=g[g.dist_cen>2.0]; k2=max(5,int(0.2*len(g2))); near2=g2.nsmallest(k2,"dist_cen"); far2=g2.drop(near2.index)
    T2=np.median(near2.pi)/np.median(far2.pi) if len(far2) and np.median(far2.pi)>0 else np.nan
    cls="acrocentric" if c in ACRO else ("metacentric" if c in META else "polycentric")
    res.append(dict(chr=c,cls=cls,n=len(g),cen_point=int(cpt),L=int(L),rho=rho,p_shift=p_shift,rho_end=rho_end,T=Tr,T_excl2Mb=T2,median_pi=float(np.median(pi))))
r=pd.DataFrame(res); r.to_csv(f"{W}/h3_perchr_{mode}.csv",index=False)
acro=r[r.cls=="acrocentric"]; poly=r[r.cls=="polycentric"]
n_sig=int(((acro.rho>0)&(acro.p_shift<0.05)).sum())
# exact permutation of class labels among the 12 non-chr10 chromosomes
vals=pd.concat([acro,poly]); T=vals["T"].values; lab=(vals.cls=="acrocentric").values
obs=T[lab].mean()-T[~lab].mean(); cnt=0; tot=0
for comb in combinations(range(len(T)),int(lab.sum())):
    m=np.zeros(len(T),bool); m[list(comb)]=True; diff=T[m].mean()-T[~m].mean(); tot+=1; cnt+= diff<=obs
p_perm=cnt/tot
# reverse direction (polycentric trough deeper than acrocentric): exact one-sided p for diff >= observed
cnt_rev=0
for comb in combinations(range(len(T)),int(lab.sum())):
    m=np.zeros(len(T),bool); m[list(comb)]=True; cnt_rev+= (T[m].mean()-T[~m].mean())>=obs
p_perm_reverse=cnt_rev/tot
T2v=vals["T_excl2Mb"].values; obs2=np.nanmean(T2v[lab])-np.nanmean(T2v[~lab])
out=dict(mode=mode,n_genes=int(len(d)),acro_sig_gradient=n_sig,mean_T_acro=float(acro["T"].mean()),mean_T_poly=float(poly["T"].mean()),perm_p=p_perm,perm_total=tot,
         mean_rho_end_poly=float(poly.rho_end.mean()),mean_rho_end_acro=float(acro.rho_end.mean()),mean_rho_acro=float(acro.rho.mean()),mean_rho_poly=float(poly.rho.mean()),chr10_T=float(r[r.cls=="metacentric"]["T"].iloc[0]) if (r.cls=="metacentric").any() else None,
         perm_p_reverse=p_perm_reverse,mean_T_excl2Mb_acro=float(np.nanmean(acro["T_excl2Mb"])),mean_T_excl2Mb_poly=float(np.nanmean(poly["T_excl2Mb"])),diff_T_excl2Mb=float(obs2))
supported = n_sig>=4 and acro["T"].mean()<0.8 and p_perm<0.05
falsified = n_sig<=2 or p_perm>0.2
telomere_confound = poly.rho_end.mean()>0 and abs(poly.rho_end.mean())>=0.5*abs(acro.rho.mean())
out["verdict"]=("supported" if supported else ("falsified" if falsified else "inconclusive"))+(" [telomere-effect confound: not claimed]" if supported and telomere_confound else "")
json.dump(out,open(f"{W}/h3_result_{mode}.json","w"),indent=1)
print(r.to_string(index=False)); print(json.dumps(out,indent=1))
