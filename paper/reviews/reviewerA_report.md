# Reviewer A — independent reproduction of Claim 1 (U. dioica tetraploid origin)

Date: 2026-10-24. Reviewer had access only to the claim, the preregistration, the code and result tables listed in the task; the author's reasoning/notebook was not read.
All numbers below were produced by code in this review (reviewerA_tier1.py, reviewerA_tier2.py) or read from the named documents.

## Verdict: **reproduced** (with one provenance gap and two methodological caveats that do not change the preregistered verdict)

## Tier 1 — statistics recomputed from gene_distances.csv without h1_test.py

Filters applied (as stated in claim/preregistration): coverage ≥ 0.8 and second-locus ratio < 0.8 in FH1, FH2, MH1, MH2 and the test target; gene on an MH2 chromosome; b and t non-missing. Parental assignment re-derived from Hirabayashi 2026 Table 2 (maternal = MH1 on chr 1,2,4,7–13; MH2 on chr 3,5,6) and checked against the d_PAT_* columns: max |Δ| = 0.

| statistic | claimed (DP) | reviewer (DP) | claimed (DA) | reviewer (DA) |
|---|---|---|---|---|
| genes | 13,439 | **13,440** | 13,420 | **13,424** |
| median b | 0.451 % | 0.449 % | 0.4566 % (JSON) | 0.455 % |
| median t | 0.497 % | 0.497 % | 0.4992 % (JSON) | 0.499 % |
| Q99(b) | 2.70 % | 2.69 % | 2.648 % (JSON) | 2.65 % |
| R = median t / median b | 1.100 (CI 1.07–1.14) | **1.106** (5-Mb block CI 1.072–1.153; gene CI 1.081–1.131) | 1.093 (CI 1.064–1.127) | **1.099** (block CI 1.067–1.144; gene CI 1.076–1.121) |
| f99 | 0.0167 (CI 0.0143–0.0199) | **0.0167** (block CI 0.0145–0.0201) | 0.0189 | **0.0190** (block CI 0.0159–0.0226) |

Bootstrap of my own design: 2,000 resamples of 5-Mb gene blocks along MH2 (111 blocks for DP, 106 for DA), plus a plain gene-level bootstrap (seed 7). Both CIs exclude R = 1 and contain the claimed point estimates.

Agreement: R within 0.006 (DP) and 0.005 (DA) of the claimed values; f99 identical to 4 decimals (DP) and 3 decimals (DA); gene counts differ by 1 (DP) and 4 (DA). Preregistered decision rule (0.67 ≤ R ≤ 1.5 and f99 ≤ 0.05) is met in every variant below.

### Robustness variants (reviewer code)

| variant | DP: n, R, f99 | DA: n, R, f99 |
|---|---|---|
| (a) max instead of min over FH1/FH2 | 13,440; R = 1.053 (1.036–1.074); f99 = 0.0146 | 13,424; R = 1.052 (1.034–1.070); f99 = 0.0161 |
| (b) t' = d(T, PAT) vs b' = d(PAT, FH1) | 13,440; R = 1.070 (1.055–1.101); f99 = 0.0135 | 13,424; R = 1.076 (1.057–1.108); f99 = 0.0153 |
| (c) ≥ 1,000 aligned columns | 7,562; R = 1.094 (1.054–1.133); f99 = 0.0159 | 7,562; R = 1.090 (1.048–1.122); f99 = 0.0179 |
| (d) chromosome 8 excluded | 12,405; R = 1.096 (1.065–1.131); f99 = 0.0164 | 12,259; R = 1.087 (1.060–1.125); f99 = 0.0190 |
| (e, added) chr 3,5,6 only (PAT = MH1; symmetric anchored-vs-anchored) | n = 2,569; R = 1.052 (gene CI 0.999–1.117) | n = 2,310; R = 1.099 (1.047–1.158) |
| (e') other 10 chromosomes (PAT = MH2 = query) | n = 10,873; R = 1.120 (1.097–1.149) | n = 11,116; R = 1.103 (1.080–1.128) |

All variants give R between 1.05 and 1.12 and f99 between 0.013 and 0.019 — the same-lineage verdict is robust to the aggregation function, the comparison design, alignment length, chromosome 8 and the PAT asymmetry.

### Provenance gap (Tier 1)
Running the author's own h1_test.py on the archived gene_distances.csv gives n = 13,438, R = 1.1062, f99 = 0.01667, Moran p = 0.052 — not the archived h1_result_DP.json values (n = 13,439, R = 1.1002, Moran p = 0.29). median(t) agrees exactly between CSV and JSON, median(b) does not (0.004490 vs 0.004515). The archived CSV (which includes BN columns) is therefore a later pipeline run than the one that produced the JSON; only one version of each file exists in the artifact store. The differences are immaterial to the verdict but the result tables and the JSON should be regenerated from one consistent run before publication.

## Tier 2 — from-scratch re-derivation on 150 random genes (seed 1)

Pipeline (independent code, reviewerA_tier2.py): downloaded GCA_048565285.1 (FH1), GCA_048565295.1 (FH2), GCA_964188135.1 (DP), GCA_054948195.1 (MH1) from the NCBI genomes FTP and the RefSeq CDS file; makeblastdb per assembly; blastn -task dc-megablast -evalue 1e-10 -max_target_seqs 5 -max_hsps 50 -dust no -soft_masking false; own HSP clustering (same subject + strand, ≤ 50 kb), own query-anchored alignment (greedy by bitscore, HSP accepted if < 10 % of its query span already covered), own p-distance over columns where both targets have ACGT (≥ 300 columns). MH2 was represented by the query CDS itself (no self-search mask). All 150 genes had hits in all four assemblies and passed cov ≥ 0.8 / ratio2 < 0.8 in my run. Genome files (3.58 GB) deleted afterwards.

| comparison (150 genes) | trim = 10 bp at HSP ends (author's choice) | trim = 0 |
|---|---|---|
| Pearson r, reviewer vs author: b / t | 0.9963 / 0.9979 | 0.9750 / 0.9917 |
| Spearman ρ: b / t | 0.9906 / 0.9952 | 0.9463 / 0.9749 |
| median |Δ| per gene: b / t | 0.00000 / 0.00000 | 0.00030 / 0.00029 |
| genes with |Δ| > 0.001: b / t | 5 / 5 | 38 / 18 |
| genes with identical value: b / t | 68.7 % / 63.3 % | 26.0 % / 24.7 % |
| subset R (reviewer) | **1.144** (gene-bootstrap CI 0.925–1.333) | 0.988 |
| subset R (author's values, same 150 genes) | 1.190 (CI 0.970–1.384) | — |
| subset f99 reviewer / author | 0.0133 / 0.0133 | 0.0133 / — |

Per-gene values are reproduced essentially exactly with the 10-bp trim (median |Δ| = 0, r ≥ 0.996). The subset-level R (1.144 vs 1.190) differs because a 150-gene median ratio is imprecise (CI width ≈ 0.4); both are consistent with the full-data 1.10.

### Caveat 1 — sensitivity to HSP-end trimming (asymmetry of b and t)
Without the 10-bp trim, subset median(b) rises from 0.00486 to 0.00556 while median(t) is unchanged (0.00550 vs 0.00556), so R falls to 0.99. Reason: on 10 of 13 chromosomes PAT = MH2 is the query itself, so b compares query bases against an anchored target (exon-boundary over-extension of an HSP puts intron bases against exon bases → mismatches), whereas t compares two anchored targets whose over-extended segments are homologous intron sequence and largely match each other. The artefact therefore inflates b more than t and biases R downward; the author's trim/self-mask removes it, and the symmetric chromosomes 3,5,6 (variant e) still give R ≥ 1.05. The direction of the artefact is conservative for the claim, but the dependence of R on an alignment-cleaning parameter should be stated in the paper.

### Caveat 2 — the Moran's I "not clustered" statement is not supported
The preregistered null (within-chromosome circular shifts of the high-t indicator) is degenerate for lag-1 adjacency: a circular shift preserves every adjacent pair except one per chromosome, so the null I ranges only 0.119–0.124 around the observed 0.1238 (DP); the reported p = 0.29 (and the 0.052 I obtain on the archived CSV) is meaningless. With a within-chromosome label-permutation null (2,000 permutations, my code), the high-t genes are strongly clustered in both DP (31 adjacent high-high pairs vs 4.1 expected; I = 0.124 vs null 99th pct 0.028; p = 0.0005) and DA (53 vs 5.1; p = 0.0005). Runs of ≥ 3 consecutive high-t genes: DP chr7 ~36 Mb (5 genes, median t = 12.9 %), chr11 ~22–23 Mb (two runs of 4, t ≈ 4.5 %), chr12 ~5 Mb (3). These are few (224 high-t genes of 13,440) and f99 = 1.7 % ≪ 0.20, so the allotetraploid branch of the decision rule is still not triggered, but the sentence "the high-divergence genes are not spatially clustered (Moran's I lag-1 p = 0.29)" must be removed or replaced by the permutation result, and the clustered high-t runs (some at U. urens-level divergence, 9–17 %) deserve a sentence (paralog/contaminant or introgressed segments).

## Document checks
- **RefSeq CDS ↔ GCA_054948205.1**: the GCF_054948205.1 assembly report states "GenBank assembly accession: GCA_054948205.1", "RefSeq assembly and GenBank assemblies identical: yes", WGS project JBQRAF01, BioSample SAMN50702769, assembly name U.dioica_male_clone48M_v1_H2. Hirabayashi 2026 Data Records: JBQRAF000000000 = MH2. All 23,131 query CDS in query_cds.fa are byte-identical to records with the same protein_id in the downloaded cds_from_genomic file, and all 48 sequence accessions in that file are in the GCF_054948205.1 report. **Confirmed.**
- **Hirabayashi 2026 Table 2** (10.1038/s41597-026-06960-7, PMC text): maternal/paternal columns read MH1/MH2 for chr 1,2,4,7,8,9,10,11,12,13 and MH2/MH1 for chr 3,5,6 — identical to the MATERNAL dict in gene_divergence.py/h1_distances.py and to the d_PAT_* columns. Male U48 is offspring of female 11-4 (Romania) × male 7-5 (Bulgaria), so PAT is the Bulgarian haplotype as claimed. **Used correctly.** Note Table 2 shows recombination on chr 4, 5, 9, 10, 12, 13 (maternal haplotype is an FH1/FH2 mosaic); this affects only the technical-floor control (min over FH1/FH2 handles it), not b.
- **DToL ploidy** (10.12688/wellcomeopenres.23399.1, Crossref TDM XML): the note states the 13 pseudomolecules represent "one subgenome of the tetraploid genome", that curation was performed "initially on all four haplotypes and subsequently on a pseudo-primary haplotype 1", BUSCO duplication 15.0 %, k-mer completeness 55 % primary / 68 % combined; the reviewer comment calls it "the first genome of a tetraploid cytotype". Flow-cytometry genome size was measured but the value is not in the text I retrieved. **Tetraploid: confirmed by the authors' statements.** The sex of specimen drUrtDioi1 (KDTOL10015, Kew, 2020-08-05) is **not stated** in the note — the claim's "tetraploid male" could not be verified from this document.

## Discrepancies summary
1. Archived gene_distances.csv does not regenerate the archived h1_result_*.json (n 13,438 vs 13,439; R 1.106 vs 1.100) — different pipeline runs. Minor, but fix before publication.
2. Moran's I test as preregistered/implemented is degenerate; high-t genes ARE clustered (permutation p = 0.0005). Claim text needs correction; verdict unchanged.
3. R depends on HSP-end cleaning (trim 0 → R ≈ 0.99 on the subset) because b and t are asymmetric with respect to alignment artefacts; the author's setting is the conservative one, but this should be disclosed.
4. Sex of the DToL specimen not verifiable from the genome note.

## Files
reviewerA_tier1.py, tier1_results.json, moran_results.json, clustering_permutation.json, symmetry_check.json, reviewerA_tier2.py, tier2_results.json, subset150.fa, subset150_author.csv, subset150_reviewer_trim10.csv, subset150_reviewer_trim0.csv.
