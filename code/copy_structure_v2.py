"""copy_structure_v2.py — pooled copy-structure baselines and Bonn contig-length context (added after independent review).
Counts, per RefSeq query, loci with CDS coverage >= 0.8 and mean HSP identity >= 97 % in each target (same rule as copy_structure.py),
then pools FH1+FH2 and MH1+MH2 (unphased-diploid baselines) and DP+DA (both tetraploid pseudo-haplotypes), and reports the
contig lengths carrying third-plus loci in the Bonn assembly. Output: results/tables/copy_structure_v2.json.
Run with the project environment: python code/copy_structure_v2.py  (≈ 10 min; streams the blast_<T>.tsv tables)."""
import sys, json, os, numpy as np
sys.path.insert(0, "code"); import gene_divergence as gd
W = "results/divergence"; os.makedirs("results/tables", exist_ok=True)
def good_loci(tg):
    out = {}; slen = {}
    for q, recs in gd.iter_query_groups(f"{W}/blast_{tg}.tsv"):
        hs = gd.parse_recs(recs)
        for h in hs: slen[h["sseqid"]] = h["slen"]
        out[q] = [l for l in gd.cluster_loci(hs) if l["cov"] >= 0.8 and l["pident"] >= 97]
    return out, slen
def dist(counts):
    c = np.array(counts)
    return dict(n=int(len(c)), frac_1=float((c == 1).mean()), frac_2=float((c == 2).mean()), frac_3=float((c == 3).mean()),
                frac_4plus=float((c >= 4).mean()), frac_ge3=float((c >= 3).mean()))
L = {}; S = {}
for t in ["FH1", "FH2", "MH1", "MH2", "BN", "DP", "DA"]:
    L[t], S[t] = good_loci(t)
res = {}
for name, a, b in [("FH1+FH2 (unphased diploid, female)", "FH1", "FH2"), ("MH1+MH2 (unphased diploid, male)", "MH1", "MH2"),
                   ("DP+DA (tetraploid Kew, both pseudo-haplotypes)", "DP", "DA")]:
    qs = set(L[a]) | set(L[b]); res[name] = dist([len(L[a].get(q, [])) + len(L[b].get(q, [])) for q in qs])
res["BN"] = dist([len(v) for v in L["BN"].values()])
third = []; best = []
for q, loci in L["BN"].items():
    if loci: best.append(S["BN"][loci[0]["sseqid"]])
    for l in loci[2:]: third.append(S["BN"][l["sseqid"]])
third = np.array(third); best = np.array(best)
res["BN_contig_len_of_3rd_plus_loci"] = dict(n=int(len(third)), median_Mb=float(np.median(third) / 1e6), frac_lt_1Mb=float((third < 1e6).mean()), frac_lt_100kb=float((third < 1e5).mean()))
res["BN_contig_len_of_best_loci"] = dict(n=int(len(best)), median_Mb=float(np.median(best) / 1e6), frac_lt_1Mb=float((best < 1e6).mean()))
json.dump(res, open("results/tables/copy_structure_v2.json", "w"), indent=1); print(json.dumps(res, indent=1))
