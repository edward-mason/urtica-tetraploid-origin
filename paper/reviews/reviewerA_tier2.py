
"""reviewerA_tier2.py — from-scratch re-derivation of per-gene p-distances for 150 random genes (Reviewer A).
Own BLAST invocation (dc-megablast, evalue 1e-10) and own anchored-alignment / p-distance code; no code shared with the author.
"""
import os, sys, gzip, subprocess, collections, json
import numpy as np, pandas as pd
BIN=r"t2\bin\ncbi-blast-2.17.0+\bin"; G="t2/genomes"; Q="t2/subset150.fa"
TG=["FH1","FH2","DP","MH1"]
def gunzip(src,dst):
    if os.path.exists(dst): return
    with gzip.open(src,"rb") as fi, open(dst,"wb") as fo:
        while True:
            b=fi.read(1<<25)
            if not b: break
            fo.write(b)
for t in TG:
    fa=f"{G}/{t}.fna"; gunzip(f"{G}/{t}.fna.gz",fa)
    if not os.path.exists(f"{G}/{t}_db.nsq") and not os.path.exists(f"{G}/{t}_db.00.nsq"):
        r=subprocess.run([BIN+r"\makeblastdb.exe","-in",fa,"-dbtype","nucl","-out",f"{G}/{t}_db","-max_file_sz","4GB"],capture_output=True,text=True)
        if r.returncode: raise SystemExit(r.stderr[-1000:])
    print(t,"db ok",flush=True)
FMT="6 qseqid sseqid bitscore qstart qend sstart send qlen qseq sseq"
hits={}
for t in TG:
    out=f"t2/blast_{t}.tsv"
    if not os.path.exists(out):
        r=subprocess.run([BIN+r"\blastn.exe","-task","dc-megablast","-query",Q,"-db",f"{G}/{t}_db","-out",out,"-outfmt",FMT,
                          "-evalue","1e-10","-num_threads","6","-max_target_seqs","5","-max_hsps","50","-dust","no","-soft_masking","false"],capture_output=True,text=True)
        if r.returncode: raise SystemExit(r.stderr[-1000:])
    rows=[l.rstrip("\n").split("\t") for l in open(out)]
    hits[t]=collections.defaultdict(list)
    for q,s,bs,qs,qe,ss,se,ql,qseq,sseq in rows:
        hits[t][q].append(dict(s=s,bs=float(bs),qs=int(qs),qe=int(qe),smin=min(int(ss),int(se)),smax=max(int(ss),int(se)),strand=int(se)>=int(ss),ql=int(ql),qseq=qseq,sseq=sseq.upper()))
    print(t,"blast ok",len(hits[t]),"queries with hits",flush=True)
def best_locus(hs,gap=50000):
    """cluster HSPs by subject+strand within `gap`; return (hsps of best-scoring cluster, score of best, score of second)"""
    groups=collections.defaultdict(list)
    for h in hs: groups[(h["s"],h["strand"])].append(h)
    clusters=[]
    for k,g in groups.items():
        g.sort(key=lambda h:h["smin"]); cur=[g[0]]
        for h in g[1:]:
            if h["smin"]-max(x["smax"] for x in cur)>gap: clusters.append(cur); cur=[h]
            else: cur.append(h)
        clusters.append(cur)
    scored=[]
    for c in clusters:
        # non-redundant score: greedy by bitscore, accept HSP if <10% of its query span already covered
        ql=c[0]["ql"]; cov=np.zeros(ql+1,bool); keep=[]; sc=0.0
        for h in sorted(c,key=lambda h:-h["bs"]):
            if cov[h["qs"]:h["qe"]+1].mean()<0.1: keep.append(h); cov[h["qs"]:h["qe"]+1]=True; sc+=h["bs"]
        scored.append((sc,keep,cov[1:].mean()))
    scored.sort(key=lambda x:-x[0])
    s2=scored[1][0] if len(scored)>1 else 0.0
    return scored[0][1],scored[0][0],s2,scored[0][2]
def anchor(keep,ql,trim):
    arr=np.full(ql,"N",dtype="<U1")
    for h in keep:
        cols=[]; p=h["qs"]-1
        for a,b in zip(h["qseq"],h["sseq"]):
            if a=="-": continue
            cols.append((p,b)); p+=1
        lo=0 if h["qs"]==1 else trim; hi=len(cols) if h["qe"]==ql else len(cols)-trim
        for p,b in cols[lo:hi]:
            if arr[p]=="N": arr[p]=b
    return arr
def pdist(a,b):
    m=np.isin(a,list("ACGT"))&np.isin(b,list("ACGT")); n=int(m.sum())
    return (int((a[m]!=b[m]).sum())/n if n>=300 else np.nan), n
qseqs={}; nm=None
for l in open(Q):
    if l.startswith(">"): nm=l[1:].strip(); qseqs[nm]=""
    else: qseqs[nm]+=l.strip()
MAT={c:("MH2" if c in (3,5,6) else "MH1") for c in range(1,14)}
auth=pd.read_csv("t2/subset150_author.csv").set_index("qid")
results={}
for trim in (0,10):
    rows=[]
    for q,qs in qseqs.items():
        ql=len(qs); anc={"MH2":np.array(list(qs),dtype="<U1")}; info={}
        for t in TG:
            if q not in hits[t]: continue
            keep,sc,s2,cov=best_locus(hits[t][q]); anc[t]=anchor(keep,ql,trim); info[t]=(cov,s2/sc if sc else np.nan,keep[0]["s"])
        chrom=int(auth.loc[q,"mh2_chr"]); pat="MH1" if MAT[chrom]=="MH2" else "MH2"
        rec=dict(qid=q,chr=chrom,PAT=pat)
        for t in TG: rec[f"{t}_cov"]=info.get(t,(np.nan,))[0]; rec[f"{t}_ratio2"]=info.get(t,(np.nan,np.nan))[1]
        def d(a,b):
            if a not in anc or b not in anc: return np.nan,0
            return pdist(anc[a],anc[b])
        for a,b in [("PAT","FH1"),("PAT","FH2"),("DP","FH1"),("DP","FH2"),("FH1","FH2"),("DP","MH1"),("DP","MH2"),("MH1","FH1"),("MH1","FH2"),("MH2","FH1"),("MH2","FH2")]:
            aa=pat if a=="PAT" else a
            rec[f"d_{a}_{b}"],rec[f"n_{a}_{b}"]=d(aa,b)
        rows.append(rec)
    df=pd.DataFrame(rows).set_index("qid"); df.to_csv(f"t2/subset150_reviewer_trim{trim}.csv")
    b=np.fmin(df.d_PAT_FH1,df.d_PAT_FH2); t_=np.fmin(df.d_DP_FH1,df.d_DP_FH2)
    ba=np.fmin(auth.d_PAT_FH1,auth.d_PAT_FH2).reindex(df.index); ta=np.fmin(auth.d_DP_FH1,auth.d_DP_FH2).reindex(df.index)
    m=b.notna()&t_.notna()
    out=dict(n=int(m.sum()),median_b=float(np.median(b[m])),median_t=float(np.median(t_[m])),R=float(np.median(t_[m])/np.median(b[m])),
             f99=float((t_[m]>np.quantile(b[m],0.99)).mean()),author_subset_R=float(np.median(ta[m])/np.median(ba[m])),author_subset_f99=float((ta[m]>np.quantile(ba[m],0.99)).mean()))
    for nm_,mine,theirs in [("b",b,ba),("t",t_,ta),("d_DP_FH1",df.d_DP_FH1,auth.d_DP_FH1.reindex(df.index)),("d_PAT_FH1",df.d_PAT_FH1,auth.d_PAT_FH1.reindex(df.index)),("d_FH1_FH2",df.d_FH1_FH2,auth.d_FH1_FH2.reindex(df.index))]:
        mm=mine.notna()&theirs.notna()
        out[f"{nm_}_pearson"]=float(np.corrcoef(mine[mm],theirs[mm])[0,1]); out[f"{nm_}_spearman"]=float(pd.Series(mine[mm].values).corr(pd.Series(theirs[mm].values),method="spearman"))
        out[f"{nm_}_medAbsDiff"]=float(np.median(np.abs(mine[mm]-theirs[mm]))); out[f"{nm_}_n"]=int(mm.sum())
    results[f"trim{trim}"]=out
    print(json.dumps(out,indent=1),flush=True)
json.dump(results,open("t2/tier2_results.json","w"),indent=1)
