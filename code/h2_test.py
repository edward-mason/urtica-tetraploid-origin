# -*- coding: utf-8 -*-

"""h2_test.py — preregistered H2 decision rule: direction of the U. dioica chr8/chr13 vs U. urens chr2 karyotype change.
Needs results/synteny/blocks_UDI_{UUR,PJU,CAN}.csv. Negative-control pair drawn with seed 20261023 from chromosomes 1-7, 9-12.
"""
import json, numpy as np, pandas as pd
S="results/synteny"; rng=np.random.default_rng(20261023)
uu=pd.read_csv(f"{S}/blocks_UDI_UUR.csv"); uu=uu[uu.block>=0]
# S8 / S13: U. dioica chr8 / chr13 genes in collinear blocks with U. urens chr2
def seg(c): return set(uu[(uu.chr_u==c)&(uu.chr_o==2)].udi)
S8,S13=seg(8),seg(13)
out={"n_S8":len(S8),"n_S13":len(S13)}
# which U. urens chromosomes do chr8 and chr13 blocks hit overall?
out["udi8_to_uur"]=uu[uu.chr_u==8].groupby("chr_o").size().sort_values(ascending=False).head(4).to_dict()
out["udi13_to_uur"]=uu[uu.chr_u==13].groupby("chr_o").size().sort_values(ascending=False).head(4).to_dict()
ctrl=sorted(rng.choice([1,2,3,4,5,6,7,9,10,11,12],2,replace=False).tolist()); out["control_pair"]=ctrl
res={}
for sp in ["PJU","CAN"]:
    b=pd.read_csv(f"{S}/blocks_UDI_{sp}.csv"); b=b[b.block>=0]
    def dist(genes):
        x=b[b.udi.isin(genes)]; return x.groupby("chr_o").size().sort_values(ascending=False)
    d8,d13=dist(S8),dist(S13)
    top8=d8.index[0] if len(d8) else None; top13=d13.index[0] if len(d13) else None
    same_chr=None; frac8=frac13=None
    if len(d8) and len(d13):
        # fraction of S8 and S13 block-genes on the single shared outgroup chromosome that maximises the minimum fraction
        cands=set(d8.index)|set(d13.index); best=None
        for c in cands:
            f8=d8.get(c,0)/d8.sum(); f13=d13.get(c,0)/d13.sum()
            if best is None or min(f8,f13)>best[0]: best=(min(f8,f13),c,f8,f13)
        same_chr=best[1]; frac8=best[2]; frac13=best[3]
    # control pair: do the two control chromosomes' blocks map to the same outgroup chromosome?
    cg=[set(uu[(uu.chr_u==c)].udi) for c in ctrl]  # genes in U. urens blocks for control chromosomes
    dc=[dist(g) for g in cg]
    ctrl_same=None
    if all(len(x) for x in dc):
        cands=set(dc[0].index)|set(dc[1].index); bestc=max(cands,key=lambda c: min(dc[0].get(c,0)/dc[0].sum(),dc[1].get(c,0)/dc[1].sum()))
        ctrl_same=(bestc,float(dc[0].get(bestc,0)/dc[0].sum()),float(dc[1].get(bestc,0)/dc[1].sum()))
    res[sp]={"S8_dist":{int(k):int(v) for k,v in d8.head(5).items()},"S13_dist":{int(k):int(v) for k,v in d13.head(5).items()},
             "shared_chr":int(same_chr) if same_chr is not None else None,"frac8_on_shared":frac8,"frac13_on_shared":frac13,
             "fission_supported":bool(frac8 is not None and frac8>=0.7 and frac13>=0.7),
             "fusion_supported":bool(frac8 is not None and (d8.iloc[0]/d8.sum()>=0.7) and (d13.iloc[0]/d13.sum()>=0.7) and top8!=top13),
             "top8":int(top8) if top8 is not None else None,"top13":int(top13) if top13 is not None else None,
             "control_pair_best_shared_chr":ctrl_same}
out["outgroups"]=res
nf=sum(r["fission_supported"] for r in res.values()); nu=sum(r["fusion_supported"] for r in res.values())
within=any(res[s]["fission_supported"] for s in ["PJU"]); within_u=any(res[s]["fusion_supported"] for s in ["PJU"])
out["verdict"]=("fission in U. dioica lineage" if (nf>=2 and within) else ("fusion in U. urens lineage" if (nu>=2 and within_u) else "inconclusive (only %d/%d outgroups agree; rule needs >=2 incl. Parietaria)"%(max(nf,nu),len(res))))
json.dump(out,open("results/h2_result.json","w"),indent=1,default=str)
print(json.dumps(out,indent=1,default=str))
