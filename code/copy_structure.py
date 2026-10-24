
import pandas as pd, numpy as np, json, sys
sys.path.insert(0,"code")
from nettle_lib import *
import gene_divergence as gd
W="results/divergence"
out={}
for tg in ["BN","DP","DA","FH1","FH2","MH1","MH2"]:
    # count near-full-length high-identity loci per query directly from BLAST output
    counts=[]; ids2=[]; n=0
    for q,recs in gd.iter_query_groups(f"{W}/blast_{tg}.tsv"):
        hs=gd.parse_recs(recs); loci=gd.cluster_loci(hs)
        good=[l for l in loci if l["cov"]>=0.8 and l["pident"]>=97]
        counts.append(len(good))
        if len(good)>=2: ids2.append(good[1]["pident"])
        n+=1
    c=np.array(counts)
    out[tg]={"n":int(n),"frac_0":float((c==0).mean()),"frac_1":float((c==1).mean()),"frac_2":float((c==2).mean()),"frac_3":float((c==3).mean()),"frac_4plus":float((c>=4).mean()),
             "second_locus_pident_median":float(np.median(ids2)) if ids2 else None,"second_locus_pident_q10":float(np.quantile(ids2,.1)) if ids2 else None}
print(json.dumps(out,indent=1))
json.dump(out,open(f"{W}/copy_structure.json","w"),indent=1)
