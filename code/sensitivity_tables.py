"""sensitivity_tables.py — small post-review tables computed from existing result files.
(1) results/tables/second_locus_identity_classes.json: per assembly, fraction of RefSeq queries with a second locus (coverage >= 0.8)
    at 80-90 %, 90-97 % and >= 97 % identity (from results/divergence/loci_<T>.csv written by gene_divergence.py summarise).
(2) results/tables/h1_residual_t_gt5_sensitivity.json: genes with t > 5 % remaining after the paralog screen, and R/f99 without them
    (from results/divergence/h1v2_genes_DP.csv / h1v2_genes_DA.csv written by h1_test_v2.py)."""
import json, numpy as np, pandas as pd, os
os.makedirs("results/tables", exist_ok=True)
out = {}
for code in ["DP", "DA", "FH1", "FH2", "MH1", "MH2", "BN"]:
    L = pd.read_csv(f"results/divergence/loci_{code}.csv"); nq = len(L); s = L[L["cov2"] >= 0.8]
    out[code] = dict(n_queries_with_best_locus=int(nq),
                     frac_2nd_cov80_pid80_90=float(((s.pident2 >= 80) & (s.pident2 < 90)).sum() / nq),
                     frac_2nd_cov80_pid90_97=float(((s.pident2 >= 90) & (s.pident2 < 97)).sum() / nq),
                     frac_2nd_cov80_pid_ge97=float((s.pident2 >= 97).sum() / nq))
json.dump(out, open("results/tables/second_locus_identity_classes.json", "w"), indent=1)
res = {}
for code in ["DP", "DA"]:
    gg = pd.read_csv(f"results/divergence/h1v2_genes_{code}.csv"); gg = gg[~gg.paralog_flag.astype(bool)]
    s = gg[gg.t <= 0.05]; hi5 = gg[gg.t > 0.05]
    res[code] = dict(n_t_gt5=int(len(hi5)), frac_t_gt5=float(len(hi5) / len(gg)), median_b_t_gt5=float(hi5.b.median()),
                     median_W_t_gt5=float(hi5.W.median()), n_excl=int(len(s)), R_excl=float(s.t.median() / s.b.median()),
                     f99_excl=float((s.t > np.percentile(s.b, 99)).mean()))
json.dump(res, open("results/tables/h1_residual_t_gt5_sensitivity.json", "w"), indent=1)
print(json.dumps(out, indent=1)); print(json.dumps(res, indent=1))
