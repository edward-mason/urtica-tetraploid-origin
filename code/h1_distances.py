
"""h1_distances.py — build per-gene pairwise p-distances from anchored sequences (preregistered common pipeline step 4-5).
Usage: python code/h1_distances.py <targets...>   (first run after summarise)
Writes results/divergence/gene_distances.csv
"""
import sys, os, json, gzip, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from nettle_lib import pdist
W=os.environ.get("NETTLE_OUTDIR","results/divergence")
MATERNAL={1:"MH1",2:"MH1",3:"MH2",4:"MH1",5:"MH2",6:"MH2",7:"MH1",8:"MH1",9:"MH1",10:"MH1",11:"MH1",12:"MH1",13:"MH1"}
targets=sys.argv[1:]
qi=pd.read_csv("results/divergence/query_index.csv").set_index("qid")
loci={t:pd.read_csv(f"{W}/loci_{t}.csv").set_index("qid").to_dict("index") for t in targets}
anch={}
for t in targets:
    with gzip.open(f"{W}/anchored_{t}.json.gz","rt") as f: anch[t]=json.load(f)
# self-derived mask: positions where the MH2 self-search anchored base equals the query base (ACGT)
def read_fa(p):
    name=None; seq=[]
    for l in open(p):
        if l.startswith(">"):
            if name: yield name,"".join(seq)
            name=l[1:].split()[0]; seq=[]
        else: seq.append(l.strip())
    if name: yield name,"".join(seq)
qseq=dict(read_fa("results/divergence/query_cds.fa"))
rows=[]
for qid,q in qi.iterrows():
    if qid not in anch.get("MH2",{}): continue
    selfa=np.array(list(anch["MH2"][qid])); Q=np.array(list(qseq[qid]))
    mask=(selfa==Q)&np.isin(Q,list("ACGT"))
    rec={"qid":qid,"gene":q.gene,"mh2_chr":q.mh2_chr,"mh2_start":q.mh2_start,"cds_len":q.cds_len,"mask_cols":int(mask.sum())}
    seqs={}
    for t in targets:
        if qid in loci[t]:
            L=loci[t][qid]
            rec[f"{t}_chr"]=L["chr"]; rec[f"{t}_seqid"]=L["sseqid"]; rec[f"{t}_start"]=L["smin"]; rec[f"{t}_cov"]=L["cov"]; rec[f"{t}_ratio2"]=L["ratio2"]; rec[f"{t}_pident2"]=L["pident2"]; rec[f"{t}_cov2"]=L["cov2"]
            arr=np.array(list(anch[t][qid]),dtype="<U1"); arr[~mask]="N"; seqs[t]=arr
    for a,b in itertools.combinations(sorted(seqs),2):
        mm,n,ind=pdist(seqs[a],seqs[b])
        rec[f"d_{a}_{b}"]=mm/n if n>=300 else np.nan; rec[f"n_{a}_{b}"]=n; rec[f"i_{a}_{b}"]=ind
    rows.append(rec)
df=pd.DataFrame(rows)
# parental labels (MH2 is the query; MATERNAL from Table 2)
chrs=df.mh2_chr
df["MAT_hap"]=[MATERNAL.get(int(c)) if pd.notna(c) else None for c in chrs]
df["PAT_hap"]=[("MH2" if m=="MH1" else "MH1") if m else None for m in df.MAT_hap]
def pick(row,a,b):
    k=f"d_{min(a,b)}_{max(a,b)}" if a!=b else None
    return row.get(k,np.nan) if k else 0.0
def dget(row,a,b):
    if not isinstance(a,str) or not isinstance(b,str): return np.nan
    if a==b: return 0.0
    k=f"d_{min(a,b)}_{max(a,b)}"; return row.get(k,np.nan)
for nm,fn in {"d_PAT_FH1":lambda r:dget(r,r.PAT_hap,"FH1"),"d_PAT_FH2":lambda r:dget(r,r.PAT_hap,"FH2"),
              "d_MAT_FH1":lambda r:dget(r,r.MAT_hap,"FH1"),"d_MAT_FH2":lambda r:dget(r,r.MAT_hap,"FH2")}.items():
    df[nm]=df.apply(fn,axis=1)
df.to_csv(f"{W}/gene_distances.csv",index=False)
print("genes",len(df)); print(df.filter(like="d_").describe().T[["count","50%","mean"]].to_string())
