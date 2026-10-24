# CLAIM LEDGER — every factual sentence of the Abstract and the Claims, with its source
Sources: J = result JSON/CSV in results/ (saved as artifacts), P = published paper (DOI), NB = notes/notebook.md entry, R = reviewer report (paper/reviews/). All numbers are copied from the named files by code (results/tables/manuscript_numbers.json, h1v2_result_*.json); nothing is from memory.

| # | Sentence (abridged) | Supporting file / field / citation |
|---|---|---|
| A1 | Diploid 2n = 26 and dominant tetraploid 2n = 52; 87 % tetraploid; uniform Cx | P Rejlová 2019 doi:10.1371/journal.pone.0218389 |
| A2 | Allopolyploid supposed (Geltman 1990); Hyb-Seq found tetraploids and diploids inseparable; no genome-wide test | P Rejlová 2019 (quote + ref 117); P Rejlová 2021 doi:10.1093/botlinnean/boaa076 (abstract); notes/literature.md; R reviewerB_report.md point 1; final novelty search NB 2026-10-24 |
| A3 | Three assemblies; parent–offspring pair; DToL = "one subgenome", HiCanu+purge_dups+YaHS, k-mer completeness 55/49/68 % | P Hirabayashi 2025 doi:10.3390/plants14010124; Hirabayashi 2026 doi:10.1038/s41597-026-06960-7 Table 2; Christenhusz 2025 doi:10.12688/wellcomeopenres.23399.1 (Table 2, Methods); ENA BioSample SAMEA7521824 (sex MALE) |
| A4 | Preregistered before download | notes/preregistration.md v1 2026-10-23T08:34:58+01:00, SHA-256 a345e28c…; data/PROVENANCE.md download log |
| A5 | n = 13,432; median t 0.50 %; median b 0.45 % | J results/divergence/h1v2_result_DP.json: paralog_filtered.n, median_t, median_b |
| A6 | R = 1.106 (chr CI 1.061–1.163; 2-Mb CI 1.077–1.142) | J h1v2_result_DP.json: paralog_filtered.R, CI_chromosomes.R_CI, CI_2Mb_blocks.R_CI |
| A7 | f99 = 1.62 % (CI 1.36 %–2.06 %); Q99(b) 2.69 % | J h1v2_result_DP.json: paralog_filtered.f99, q99_b; CI_chromosomes.f99_CI |
| A8 | Held-out DA: R = 1.097 (CI 1.065–1.137); f99 = 1.79 % (CI 1.45 %–2.20 %); n = 13,404 | J h1v2_result_DA.json |
| A9 | High-t genes clustered: DP 30 vs 3.9 (p 0.0005); DA 46 vs 4.5 (p 0.0005) | J h1v2_result_DP/DA.json: clustering |
| A10 | Blocks chr7 0.8–22.8 Mb (t 0.62 %, b 0.16 %, W 0.15 %), chr8 10.9–13.9 Mb (t 0.58 %, b 0.16 %); null99 0.34 %; max windowed median t 1.25 % | J h1v2_result_DP.json: paired_scan.blocks, null99_rolling_median_diff, max_windowed_median_t; Figure 5 |
| A11 | d(DP,DA) median 0.33 %, 2.9 % > Q99(b); log-scale dip p 0.99; W p 0.98; b p 1.00; raw p 0 | J h1v2_result_DP.json: dip |
| A12 | Second loci ≥ 97 %: 2,639/16,538 genes, identity median 98.96 %; 90–97 % class 10.3/9.9 vs 9.4–9.8 %; 80–90 % 4.9–5.2 vs 5.4–5.6 %; ≥ 97 % 22.3/21.4 vs 9.6–12.2 % | J results/tables/manuscript_numbers.json (dp_dup_*); results/tables/second_locus_identity_classes.json |
| A13 | Power model exclusions (δ ≥ 0.5 % at φ ≥ 0.25; δ ≥ 2 % at φ ≥ 0.10; δ = 0.25 % at φ ≥ 0.5) | J h1v2_result_DP.json: power_model, excluded_combos; Table 3 |
| A14 | Symmetric chr 3/5/6 R 1.052 vs others 1.119 | J h1v2_result_DP.json: symmetric_chr356, other_chr |
| A15 | Trim-0 sensitivity R 1.114, f99 1.74 %, n 13,539 | J results/tables/h1_trim0_sensitivity.json |
| A16 | Residual t > 5 % genes 80/82; excluding them R 1.100/1.093, f99 1.32 %/1.51 % | J results/tables/h1_residual_t_gt5_sensitivity.json |
| A17 | Floor exactly 0 in 98.8 % of 16,819 genes; W median 0.47 %, q99 3.42 %; b median 0.45 %, q99 2.77 %; d(FH1,UU) median 9.4 % | J manuscript_numbers.json; gene_distances.csv |
| A18 | Bonn: ≥ 3 loci 23.9 %, 1 locus 17.1 %, 2 loci 51.0 %; contigs median 40 Mb; pooled baselines 11.5 %, 13.3 %, 32.1 %; R_BN single-copy 1.108 (CI 1.03–1.18), best-locus 0.940; d(BN,DP) 0.59 % | J results/tables/copy_structure_v2.json, copy_structure.json, h1_result_BN.json, h1_BN_bestlocus.json; P Wolff 2026 doi:10.64898/2026.05.15.725508 ("two haplophases", single root, Braunschweig) |
| A19 | H3 female: 5/5 acro ρ > 0; T acro 0.67 vs poly 0.38; chr10 0.40; p 1.00 / reverse 0.001; end ρ -0.37/-0.11; excl-2Mb 0.74/0.41 | J results/tables/h3_result_female.json, h3_perchr_female.csv |
| A20 | H3 male: T 0.77 vs 0.44; reverse p 0.018; excl-2Mb 0.75/0.48; chr13 interval at opposite end | J h3_result_male.json, h3_perchr_male.csv; P Hirabayashi 2026 Supplementary Table 5; Hirabayashi 2025 Table S8 |
| A21 | H2: RBH 14,214/13,635/13,189; S8 1,252, S13 739; PJ 453/423, 263/219; CAN 449/151/147, 159/147; control pair disjoint; eight segments; 5 telomeric arrays within 11 kb of start | J results/tables/h2_result.json, h2_itr.json, blocks_UDI_*.csv |
| A22 | Ramie/Morus dropped: no annotated chromosome-level assembly on NCBI | data/annotations.csv; data/assemblies.csv |
| A23 | Reviewer A reproduced (r 0.996–0.998) ; Reviewer B points addressed | R reviewerA_report.md, reviewerB_report.md; paper/REVIEWER_RESPONSES.md |
