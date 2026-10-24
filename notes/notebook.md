# Lab notebook — Urtica dioica comparative genomics

## 2026-10-22T21:48:16+01:00 — Phase 0: workbench

- Environment: Claude Science sandbox on Windows 11 host (AMD Ryzen 7 3700U, 4C/8T, 16 GB RAM; sandbox sees 14 GiB, ~4 GiB free at start). No GPU.
- WSL 2 probe: `wsl --status` → `Wsl/E_ACCESSDENIED` from inside the sandbox. Decision (user-approved): Windows-native pipeline. Read-level analyses (alignment of RNA-seq / resequencing / methylation reads) are OUT OF SCOPE; only published processed outputs may be used.
- Disk at start: 118.3 GB free on workspace drive. Hard peak ceiling 80 GB.
- Network: NCBI Datasets/E-utilities/FTP, ENA portal, Ensembl Plants REST, CrossRef, Europe PMC reachable. Dryad (datadryad.org) blocked by allowlist — will request only if a Dryad-only file is needed. All 26 attached literature connectors disconnected; using direct OpenAlex/CrossRef/Europe PMC APIs instead.
- Conda env `nettle` (py3.11): biopython, pyhmmer, pyfastx, statsmodels, scikit-learn, pandas, numpy, scipy, matplotlib, seaborn, snakemake-minimal, parasail-python. Creation started 2026-10-22T21:48:16+01:00.
- External Windows binaries to fetch into code/bin: NCBI BLAST+, MAFFT, IQ-TREE 2, PAML.

### Disk / runtime log
| step | est. peak disk | est. runtime | actual |
|---|---|---|---|
| env creation | <3 GB | ~10 min | |

### 2026-10-22T22:27:44+01:00 — Phase 0 outcome

Tooling test results inside the Claude Science sandbox:
- NCBI BLAST+ 2.17.0 (x64-win64 tar.gz, MD5 verified against NCBI .md5): RUNS.
- PAML 4.10.10 (win-x86_64): codeml.exe starts (smoke test 2026-10-22). Correction 2026-10-23: baseml and yn00 were NOT tested on 2026-10-22; tested now — baseml: 'BASEML in paml version 4.10.10, 27 Jan 2026'; yn00: 'YN00 in paml version 4.10.10, 27 Jan 2026' (each invoked with a missing control file; a start-up message means the binary loads).
- IQ-TREE 2.4.0, 2.3.6, 1.6.12 Windows builds: all fail with STATUS_ILLEGAL_DLL_RELOCATION (WinError 623) — a sandbox loader mitigation against the MinGW-built binaries; not fixable from here. Removed from code/bin.
- MUSCLE 5.3 win64: same WinError 623. Removed.
- MAFFT: mafft.cbrc.jp outside allowlist; not requested (pyfamsa covers MSA).
- No win-64 builds exist on conda-forge for iqtree/fasttree/raxml-ng/mafft/muscle/hmmer/diamond/mmseqs2/minimap2/samtools/bedtools/mcscanx/orthofinder (checked conda-forge win-64 repodata, 2026-10-22).

Decisions:
- MSA: pyfamsa (FAMSA, Windows wheel) for protein; codon alignments by back-translation in Biopython.
- Pairwise alignment / homology search: BLAST+ (tblastn/blastp/blastn), pyhmmer (HMMER3 via Python), Biopython PairwiseAligner, pyopal.
- Tree inference: ML via PAML baseml/codeml (fixed topology from NJ or stepwise search for small trees); NJ via scikit-bio/Biopython. Large ML phylogenomics is OUT OF SCOPE; hypotheses must not depend on it.
- dN/dS: PAML codeml / yn00.
- Synteny: BLAST-based reciprocal best hits + custom collinearity chaining in Python (MCScanX unavailable).
- Repeat annotation: reuse published annotations only; RepeatModeler impossible here.
- Environment `nettle` is a pip-backed venv (python 3.11): biopython, pyhmmer, pyfastx, statsmodels, scikit-learn, pandas, numpy, scipy, matplotlib, seaborn, snakemake-minimal, parasail-python, pyfamsa, biotite, scikit-bio, pyskani, pyopal, gffutils, intervaltree, pyarrow.

Disk/runtime actuals: env build ~30 min wall (package downloads), binaries ~150 MB on disk; peak disk <1 GB.

### 2026-10-22T22:33:30+01:00 — User directives
- User granted blanket approval to run all code and to proceed autonomously overnight without asking. Network domains granted: datadryad.org, zenodo.org, figshare.com, ndownloader.figshare.com, static-content.springer.com, www.mdpi.com, download.cncb.ac.cn (ftp.ebi.ac.uk, ftp.ensemblgenomes.ebi.ac.uk, rest.uniprot.org, api.biorxiv.org, api.semanticscholar.org, doi.org already allowed). OpenAlex credential approved.
- The brief's stop conditions (80 GB peak disk, 72 h runtime, GPU, paid compute, irreversible actions) remain in force: I will log and halt at those, not ask for routine steps.
- Phase 1 dispatched as three parallel sub-agents: Genomes, Reads, Literature (2026-10-22T22:33:30+01:00).

## 2026-10-23T08:42:14+01:00 — Phase 3 start
- Preregistration saved (artifact 8f9e6445, sha256 a345e28c…) at 2026-10-23T08:34:58+01:00 before any FASTA download.
- Downloads: NCBI FTP (FH1, FH2, MH1, MH2+RefSeq annotation, DP, DA, UU, PJ; Morus/Cannabis proteins+GFF), Ensembl rapid genesets UU/PJ 2024_02, figshare 230981 (BRAKER3 gene sets female/male, EDTA female TE annotation, male supplement PDF), Plants 2025 supplement via Europe PMC PMC11722821 (MDPI returns 403 to automated clients). Figshare downloads redirect to s3-eu-west-1.amazonaws.com (allowlisted on request).
- Centromere operational rule (fixed before any divergence value was computed): Table S8 intervals < 100 kb ignored (the recurrent 1,250,000 starts are an analysis-window edge); acrocentric chromosome point = midpoint of the largest interval; polycentric/metacentric point = midpoint of the span of all remaining intervals. Saved data/annot/centromere_points.csv.

### 2026-10-23T09:18:59+01:00 — BLAST parameters
- Benchmark (400 CDS vs MH2 db, 8 threads): megablast 10 s; dc-megablast 89 s (max_target_seqs 5, max_hsps 40). First attempt with max_target_seqs 10 / max_hsps 60 projected 2.6 h per genome and was stopped at 3,875/23,131 queries.
- Decision: dc-megablast for ALL seven targets (identical method for baseline and test comparisons; megablast word-28 seeding could under-align the more divergent targets and bias t vs b asymmetrically). Estimated 90 min per genome, ~10.5 h total; outputs ~0.4 GB each. Preregistration said "task blastn"; dc-megablast is the discontiguous-seed variant of the same search — recorded as a parameter deviation, not a statistical one.

### 2026-10-23T12:20:48+01:00 — pipeline calibration on the MH2 self-search (not a hypothesis test)
- Self-search (RefSeq CDS vs its own genome): 23,131 queries all hit; best locus at the annotated position for 94.8 %; chromosome concordance 97.7 % after stripping BLAST's `ref|…|` seqid wrapper.
- Self p-distance floor with no trimming: mean 0.00079, 22 % of genes > 0.001, 99th pct 0.8 % — exon-boundary over-extension of HSPs into intron sequence. HSP end-trim calibration (1,500 genes): trim 0/5/10/15/20 → mean 0.00079/0.00054/0.00026/0.00020/0.00018; coverage 0.98/0.96/0.93/0.90/0.87. Decision: trim = 10 (fixed now, before any between-haplotype distance is seen) PLUS a per-gene mask excluding query positions where the self-search anchored base ≠ query base. All pairwise distances use only unmasked columns where both targets have ACGT.
- Summarise step rewritten without pandas indexing (0.74 s → 0.012 s per query).

## 2026-10-23T17:39:35+01:00 — H3 discovery result (female data; held-out male data NOT yet examined)
Baseline controls (16,819 genes pass 4-haplotype single-copy filter): technical floor median 0 (99.7 % of genes exactly 0); W = d(FH1,FH2) median 0.474 %, q95 1.77 %, q99 3.42 %; b = min d(PAT,FH1/FH2) median 0.452 %, q95 1.53 %, q99 2.77 %. Within-individual ≈ between-population divergence (Romania vs Bulgaria).
H3 verdict by preregistered rule: **falsified** (n = 16685 genes).
- Acrocentric chromosomes: 5/5 show rho > 0 with shift-null p < 0.05; mean trough ratio T = 0.668 (< 0.8 criterion met).
- Polycentric chromosomes: mean T = 0.377 — a DEEPER central trough than the acrocentric class; exact label permutation p = 1.000 (792 assignments); chr10 (metacentric) T = 0.399.
Interpretation (written before held-out data): CDS nucleotide diversity is depressed around the RepeatOBserver-predicted centromeric region on EVERY chromosome, most strongly on the seven "polycentric" ones, whose predicted region is a broad central span. The landscape is that expected of a broad central low-recombination domain, not of distributed recombination; it does not, however, discriminate polycentric from monocentric organisation because distally biased crossovers also occur in holocentrics (Castellani 2024). H3 is reported as a negative result. The per-chromosome table is results/divergence/h3_perchr_female.csv.
Caveat: the predicted regions coincide with gene-density minima (Table S8) and with TE-rich DNA; CDS diversity in gene-poor regions may also be reduced by background selection acting over larger linked distances.

## 2026-10-23T18:32:21+01:00 — H1 discovery result (DP = GCA_964188135.1; held-out DA NOT yet examined)
Preregistered statistics (code/h1_test.py, seed 20261023), n = 13439 genes passing single-copy/coverage filters in FH1, FH2, MH1, MH2 and DP:
- median b = 0.00451, median t = 0.00497, median W = 0.00468; R = 1.100 (2-Mb block bootstrap 95 % CI 1.073–1.136).
- Q99(b) = 0.0270; f99 = 0.0167 (CI 0.0143–0.0199); f95 = 0.0681 (CI 0.0634–0.0738).
- Moran's I (lag 1, high-t indicator) = 0.1238; null 99th pct = 0.1239; p = 0.291 → no spatial clustering.
- Technical floor ok (0.0). Scale control pending U. urens.
**Verdict by the preregistered rule: same-lineage tetraploid** (0.67 ≤ R ≤ 1.5 ✓; f99 ≤ 0.05 ✓; d(DP,DA) bimodality to be checked when DA is opened).
Artifact checks:
- Duplicated loci in DP: 2,639/16,538 genes (16.0 %) have a second DP locus ≥ 97 % identity, ≥ 80 % coverage (diploid haplotypes: 6–7 genes) — the expected signature of retained homoeologs/homologs in a tetraploid-derived pseudo-haplotype. Their best-locus divergence from the female (median 0.0050) equals that of single-copy genes (0.0050), and the second loci themselves are 98.96 % identical (median; q10 97.94 %) to the diploid CDS — no hidden divergent subgenome among duplicated loci.
- 440 genes whose best DP locus is on an unplaced scaffold; 20 chromosome mismatches; R without chr8 = 1.091, f99 = 0.0165.
- Per-chromosome R 0.99–1.33; f99 ≤ 0.03 on every chromosome. The chr7 (0–30 Mb) and chr8 (10–20 Mb) R inflation arises from LOW baseline b (0.0013 / 0.0010; 26–30 % of genes b = 0), i.e. the Bulgarian paternal haplotype is near-identical to one female haplotype there, while t is normal (0.0034 / 0.0067). Chr8 10–20 Mb overlaps the inversion-bearing region — a shared low-diversity haplotype block consistent with suppressed recombination; recorded as an observation only (chr8 is the UBC group's SDR candidate).
- Quantiles: t q99 = 0.0359 vs b q99 = 0.0270 and W q99 = 0.0322 — the tetraploid's upper tail lies within the diploid within-individual range.
- DP's nearest female haplotype coincides with PAT's nearest female haplotype in 74.4 % of genes (50 % expected if independent) — the same female haplotype tends to carry the more common allele for both the Bulgarian diploid and the UK tetraploid.
Interpretation: the DToL tetraploid pseudo-haplotype is a sample from the same gene pool as the Eurasian diploids; nothing in 13,439 genes indicates a divergent second parental lineage. Awaiting DA (held-out replicate), U. urens (scale) and BN (exploratory).

## 2026-10-23T20:00:53+01:00 — H1 held-out replicate (DA = GCA_964188125.1) and scale control
- Scale control: median d(FH1, U. urens) = 0.0930 (≈ 20 × median b) → dynamic range adequate (criterion ≥ 3×).
- Held-out DA (n = 13420): R = 1.093 (CI 1.064–1.127); f99 = 0.0189 (CI 0.0159–0.0223); f95 = 0.0679. Verdict: **same-lineage tetraploid** — identical to DP → H1 same-lineage verdict CONFIRMED by the preregistered held-out rule.
- DA Moran's I = 0.194 (p = 0.006): the 254 genes above Q99(b) show weak spatial clustering (largest cluster 13 genes in chr11 20–25 Mb), but f99 = 1.9 % is an order of magnitude below the 20 % allotetraploid criterion. For comparison the diploid's own within-individual W exceeds Q99(b) in 1.533% of the same genes.
- Secondary (within-tetraploid) d(DP,DA): n = 10,956, median 0.0031 (below W = 0.0047), fraction > Q99(b) = 2.3–2.9 %. Hartigan dip on raw p-distances rejects unimodality for d(DP,DA) AND for W AND for b alike (p ≈ 0) because of the spike at zero (discrete, zero-inflated distances); on log10 of non-zero values all three are unimodal (p = 0.99). Conclusion: no homoeologous second mode; the primary–alternate pairs behave like allelic pairs of a diploid. (The dip component of the preregistered rule is therefore evaluated on the log scale with the diploid baselines as reference; the raw-scale test is non-diagnostic and this is declared as a post-hoc clarification.)
- Pipeline note: DA has no chromosome-level scaffolds, so DA_chr is undefined for all genes; the DA replicate therefore also serves as a check that the result does not depend on chromosome assignment.

## 2026-10-23T23:24:46+01:00 — Exploratory: copy structure of the Bonn assembly (Amendment 1) and synteny sanity check
Loci per RefSeq CDS with coverage ≥ 0.8 and identity ≥ 97 % (code/copy_structure.py), all queries:

| target | n queries | 1 locus | 2 loci | 3 loci | ≥4 loci | 2nd-locus identity (median %) |
|---|---|---|---|---|---|---|
| BN | 23125 | 0.171 | 0.510 | 0.158 | 0.082 | 99.11 |
| DP | 23016 | 0.663 | 0.136 | 0.025 | 0.063 | 98.88 |
| DA | 22898 | 0.665 | 0.134 | 0.020 | 0.061 | 98.91 |
| FH1 | 23111 | 0.814 | 0.034 | 0.012 | 0.058 | 98.73 |
| FH2 | 23112 | 0.806 | 0.034 | 0.010 | 0.054 | 98.60 |
| MH1 | 23112 | 0.801 | 0.034 | 0.011 | 0.056 | 98.70 |
| MH2 | 23131 | 0.849 | 0.048 | 0.015 | 0.061 | 98.88 |

- Diploid haplotype assemblies (FH1, FH2, MH1, MH2) carry one locus for 80–85 % of genes; ≥ 3 loci for 6.6–7.7 % (gene families).
- DToL tetraploid pseudo-haplotypes (DP, DA): two loci for 13.4–13.6 % (retained homologs/homoeologs), ≥ 3 for 8.1–8.8 %.
- Bonn assembly (BN, 1.1 Gb hifiasm primary): two loci for 51.0 %, three for 15.8 %, ≥ 4 for 8.2 %; only 17.1 % single-copy. The excess of genes with ≥ 3 near-identical loci (23.9 % vs 6.6–7.7 % in diploid haplotypes, 8.1–8.8 % in DP/DA) cannot arise from a diploid (max two allelic copies) and indicates ≥ 3 resolved haplotypes at ~16–17 % of loci — the signature of a partly collapsed tetraploid. Together with the 1.1 Gb span and the 1.3 Gb mapping-based size estimate, the Bonn plant is most parsimoniously a tetraploid; the preprint's "diploid" inference rests on a two-copy BUSCO peak, which is exactly what partial collapse of four haplotypes to ~two would produce. Reported as an observation with the caveat that no flow cytometry exists for that plant.
- Synteny sanity check (U. dioica MH2 RefSeq vs U. urens Ensembl 2024_02; 14,214 RBH pairs, 13,418 in 198 collinear blocks): every U. dioica chromosome maps 1:1 to a single U. urens chromosome, except chr8 (1,252 block genes) and chr13 (744) which both map to U. urens chr2 — reproducing Hirabayashi 2025 Figure 5.

## 2026-10-24T00:10:37+01:00 — Exploratory BN divergence (Amendment 1; not part of the H1 decision)
- Preregistered statistic restricted to BN single-copy genes (n = 2706): R = 1.108 (CI 1.027–1.185), f99 = 0.0181 → would read "same-lineage".
- Best-locus (nearest BN haplotype copy) over all genes passing the diploid filters (n = 16604): R = 0.940, f99 = 0.0118; multi-copy genes only (n = 13897): R = 0.905. R < 1 is expected when the minimum is taken over 2–4 resolved haplotypes instead of one.
- median d(BN, DP) = 0.0059 (Germany tetraploid vs UK tetraploid) vs W = 0.0047: tetraploid–tetraploid divergence is at the diploid allelic level.
- Conclusion: both tetraploid-candidate individuals (Kew DToL, Bonn) sample the same gene pool as the Eurasian diploids.

## 2026-10-24T02:03:35+01:00 — H3 held-out replicate (male: pi = d(MH1,MH2), MH2 coordinates, male Supplementary Table 5 centromeres)
Verdict: **falsified** (n = 16654): acrocentric gradient in 4/5; mean T(acro) = 0.766, mean T(poly) = 0.440; exact permutation p = 0.984; chr10 T = 0.303. Identical conclusion to the female discovery set: the within-individual diversity trough around the predicted centromeric span is deeper on the "polycentric" chromosomes than on the acrocentric ones.
Chr13 anomaly: the male Table S5 interval (1.25–1.455 Mb) lies at the opposite chromosome end from the female Table S8 interval (18.9–19.6 Mb of 20.9 Mb); in the male, diversity is HIGHER near that point (T = 1.36, rho = −0.18). The 1.25 Mb start recurs across chromosomes in both tables and is almost certainly a RepeatOBserver analysis-window edge; the chr13 centromere call in the male should be treated as unresolved. Removing chr13 does not change the class comparison (polycentric trough still deeper).
Pipeline note: subprocesses launched from the stdlib control kernel need the conda env's Library\bin on PATH, otherwise numpy's BLAS DLL fails to delay-load (exit 0xC06D007F) — fixed by setting PATH explicitly.

## 2026-10-24T06:46:18+01:00 — H2 result (direction of the chr8/chr13 – U. urens chr2 change)
Synteny (RBH evalue 1e-10, ≥50 % mutual coverage; blocks ≥5 pairs within 20-gene windows): U. urens 14,214 RBH / 13,333 in 305 blocks; Parietaria 13,635 / 10,397 in 642 blocks; Cannabis 13,189 / 7,936 in 605 blocks. Every U. dioica chromosome maps 1:1 to one U. urens chromosome except chr8 (1,252 block genes) and chr13 (739), both to U. urens chr2 — reproducing Hirabayashi 2025.
**Verdict by the preregistered rule: inconclusive (only 0/2 outgroups agree; rule needs >=2 incl. Parietaria).** In Parietaria, S8 genes split PJ1 (453) / PJ3 (423) and S13 genes split PJ5 (263) / PJ6 (219); in Cannabis, S8 → CAN1 (449), CAN4 (151), CAN8 (147) and S13 → CAN10 (159), CAN3 (147). No outgroup places ≥70 % of either segment on one chromosome, so neither direction can be called. Negative-control pair (chr4, chr6; seed 20261023) maps to disjoint outgroup chromosomes in both outgroups (shared fractions ≤ 0.08), confirming the test can return "separate".
Post-hoc observations (exploratory, not part of the decision):
- In BOTH outgroups the chr8-derived and chr13-derived gene sets occupy entirely disjoint chromosomes (0 % of S13 genes on the chromosome carrying most S8 genes). Under "fission in U. dioica" (ancestral single chromosome) both outgroups would have had to split the ancestral chromosome independently; under "fusion in U. urens" no outgroup event is required. Parsimony therefore leans to fusion in the U. urens lineage, but the preregistered threshold is not met and this is reported as a lean, not a result.
- Parietaria alternates between chromosome pairs (PJ1/PJ3 along chr8, PJ5/PJ6 along chr13) at fine scale and has 35,329 genes vs 18,378 in U. urens — the signature of retained homoeologs after a lineage-specific polyploidy; it is unsuitable as a chromosome-identity outgroup. Cannabis shows deep Rosales-level rearrangement (chr8 ↔ CAN1 + CAN4 + CAN8; chr13 ↔ CAN3 + CAN10).
- U. urens chr2 (32.2 Mb) is a MOSAIC of eight alternating segments: chr8 0.27–0.51 | chr13 7.2–18.7 | chr8 45.2→35.2 (inverted) | chr13 1.1–3.5 | chr8 2.4–35.1 (831 genes, intact) | chr13 0.6–7.1 | chr8 0.5–2.4 | chr13 3.5–4.5 Mb. Seven switch points imply a fusion/fission plus ≥3 inversions spanning the junction in one lineage or the other — the relationship is more complex than a single fusion/fission event.
- Telomeric repeat arrays (TTTAGGG)n≥4 on U. urens chr2: 5 arrays, all within the first 11 kb (start telomere); no interstitial arrays — no ITR support for a simple fusion junction (which the mosaic structure would in any case have dispersed).

## 2026-10-24T08:10:51+01:00 — Phase 4: independent review and corrections
Reviewer A (reproduction; artifact reviewerA_report.md 7f6160b0): Tier 1 own code R = 1.106/1.099, f99 = 0.0167/0.0190 (DP/DA) — matches. Tier 2 from-scratch BLAST on 150 random genes: per-gene r = 0.996–0.998, median |Δ| = 0. Findings accepted: (i) archived CSV and JSON came from different pipeline runs (MH2 loci regenerated after the JSON) — FIXED by regenerating all statistics from the final gene_distances.csv (h1_test_v2.py); (ii) circular-shift Moran's I null is degenerate (shifting the whole indicator preserves adjacency) — REPLACED by within-chromosome permutation of adjacent high-t pairs; (iii) R sensitive to HSP end-trim on the 150-gene subset — full-data trim-0 sensitivity run (results/divergence_trim0, below).
Reviewer B (hostile; artifact reviewerB_report.md 02df3509): 18 points. Accepted and acted on:
- DToL pipeline is HiCanu + purge_dups + YaHS (not hifiasm) — corrected in text.
- Rejlová et al. (Bot J Linn Soc 195:606–621; volume year 2021, online 2020) already showed tetraploids coalesce with diploids by Hyb-Seq — framing changed to "first genome-wide, assembly-based test"; Geltman 1990 named as source of the allopolyploid hypothesis via Rejlová 2019.
- Spatial clustering of high-t genes IS significant: DP 30 adjacent pairs vs 3.9 expected (p = 0.0005); DA 46 vs 4.5 (p = 0.0005). Windowed paired (t − b) scan, 50-gene rolling median vs within-chromosome permutation 99th pct (0.34 pp): DP blocks chr7 0.8–22.8 Mb (t 0.62 %, b 0.16 %, W 0.15 %, 22 % b = 0); chr8 10.9–13.9 Mb (t 0.58 %, b 0.16 %, W 0.85 %, 20 % b = 0); DA blocks chr2 31.9–38.4 Mb (t 0.52 %, b 0.21 %, W 0.53 %, 26 % b = 0); chr8 12.4–13.1 Mb (t 0.50 %, b 0.15 %, W 0.92 %, 24 % b = 0). In every block b is deflated (shared low-diversity diploid haplotype), not t elevated; max 50-gene median t anywhere = 1.25 % (DP), 1.31 % (DA).
- Paralog filter (d(T,UU) − d(FH1,UU) > 5 pp): flags 6 (DP) and 15 (DA) genes; filtered R = 1.106 / 1.097, f99 = 0.0162 / 0.0179. The 82 DA genes with t > 5 % are not paralogs in the filtered set: DP shows the same high t (median 5.7 %) and all haplotypes are ~2× more divergent from U. urens than average (fast-evolving, highly polymorphic genes).
- Chromosome-level bootstrap CIs: R 1.061–1.163 (DP), 1.065–1.137 (DA); f99 0.0136–0.0206 / 0.0145–0.0220.
- Dip test now executed in-script: raw p = 0 for d(T,other), W and b alike; log-scale p = 0.99/0.98/1.00. Symmetric-design chromosomes 3/5/6: R = 1.052 vs 1.119 elsewhere.
- Power model (b resampled; fraction φ of genes + δ): combinations whose expected R or f99 exceed the observed chromosome-level upper CI are excluded: all δ ≥ 0.5 % at φ ≥ 0.25; δ ≥ 2 % at φ ≥ 0.10; δ = 0.25 % at φ ≥ 0.5. NOT excluded: δ ≤ 1.5 % at φ = 0.10, δ = 0.5 % at φ = 0.10, δ = 0.25 % at φ ≤ 0.25. Claim 1 restated as the excluded region.
- k-mer completeness of DToL (55 % primary / 68 % combined vs 98 % for the diploid pair) acknowledged as the principal limitation: an undeposited divergent haplotype pair cannot be excluded by assembly-only analysis; purge_dups places purged haplotigs in the alternate (DA), which shows no divergent class, but DA is itself 49 % complete. Read-based k-mer analysis (GenomeScope2 -p 4 / Smudgeplot on ERR10662027–9) is required and listed in VALIDATION.md; infeasible here (no WSL/Linux; pure-Python k-mer counting of 45 Gb HiFi not practical).
- Bonn baseline: unphased-diploid baselines ≥3 loci = 11.5 % (FH1+FH2), 13.3 % (MH1+MH2); tetraploid DP+DA = 32.1 %; BN = 23.9 %. BN's third-plus loci lie on contigs of median 40.3 Mb (0.06 % on contigs < 1 Mb). Claim 2 softened to "consistent with a partly collapsed tetraploid; haplotig triplication in a diploid not excluded without read depth/flow cytometry".
- H3: opposite-direction effect reported: exact one-sided p(T(poly) < T(acro)) = 0.001 (female), 0.018 (male). Preregistered controls: distance-to-nearest-end rho = -0.37 (poly) vs -0.11 (acro) in the female; after removing genes within 2 Mb of the centromere point T(acro) = 0.74 vs T(poly) = 0.41 (female), 0.75 vs 0.48 (male).
- H2: Morus notabilis (RefSeq GCF_000414095.1 is scaffold-level; chromosome-level Morus assemblies carry no NCBI annotation) and ramie (no NCBI annotation for GCA_018132145.1; CNCB GWH annotation not retrieved) were dropped — now declared as a deviation; Parietaria WGD (35,329 vs 18,378 genes) makes the ≥70 % rule unattainable, so "inconclusive" partly reflects outgroup choice; Hirabayashi 2025 already reported the chr8/13–chr2 synteny; the new element is the eight-segment structure.
Rebutted / not actionable: (B-15c) tissue pooling in Wolff et al. — Methods not checkable from the PDF text layer here; stated as unknown. (B-10) thresholds "unjustified by a model" — accepted; model now supplied, verdict reframed.


## 2026-10-24T08:59:03+01:00 — Phase 5: post-review revision, novelty search, packaging

**Full-data trim sensitivity.** Re-ran summarise/distances/h1 with NETTLE_TRIM=0 (results/divergence_trim0): n = 13,539, R = 1.114, f99 = 1.74 % vs 1.106 / 1.62 % at trim 10. Reviewer A's subset-level discrepancy (150 genes) was sampling noise. Saved results/tables/h1_trim0_sensitivity.json.

**Figure 5 (new).** 50-gene rolling medians of W, b, t, tA, d(DP,DA) along all 13 MH2 chromosomes with the two paired-scan blocks shaded. First draft had per-panel x-scales with labels only on the bottom panel (misleading) — fixed to a common 0–62 Mb axis with chromosome-end ticks. Observation recorded from the figure: the tetraploid tracks overlie the diploid tracks everywhere; in the chr7 block DP carries a 0.6 % haplotype where FH1, FH2, PAT and DA share a ~0.1 % one; d(DP,DA) spikes > 1 % in short windows on chr1 (~53 Mb), chr3 (~14, ~31 Mb), chr4 (30–36 Mb), chr5 (~17, ~23 Mb) — in each, one pseudo-haplotype is near-identical to the female and the other at ordinary allelic divergence (0.6–0.7 %); no window has both t and tA > 1.3 %.

**Residual high-divergence genes.** After the paralog screen, 80 (DP) and 82 (DA) genes keep t > 5 % (0.6 % of each set): hypervariable genes (median b 1.0 % / 1.9 %; U. urens distance 17–19 % vs 9.3 % typical), spread over all chromosomes; d(DA,UU) ≈ d(F,UU) so the > 5-point paralog rule cannot flag them. Excluding them: R 1.100/1.093, f99 1.32 %/1.51 %. Saved results/tables/h1_residual_t_gt5_sensitivity.json. Noted that DA has no chromosome assignments at all (contig-level), so "unplaced" is uninformative for DA.

**Second-locus identity classes (own recomputation of Reviewer B point 7).** 90–97 %: DP 10.3, DA 9.9, FH1 9.5, FH2 9.8, MH1 9.8, MH2 9.4 %; 80–90 %: 4.9/5.2 vs 5.4–5.6 %; ≥ 97 %: 22.3/21.4 vs 9.6–12.2 % (BN 74.8 %). Saved results/tables/second_locus_identity_classes.json.

**Manuscript rewritten** (paper/manuscript.md, ~54 kB): title no longer asserts autotetraploidy; Abstract/Results/Discussion carry h1_test_v2 numbers (n 13,432, R 1.106 chr-CI 1.061–1.163, f99 1.62 % CI 1.36 %–2.06 %; DA R 1.097, f99 1.79 %); clustering reported as positive and decomposed; power model Table 3; Bonn softened to "probable"; H3 restated with opposite-direction p and both controls; H2 deviations declared; Rejlová 2021 as prior molecular evidence; DToL pipeline corrected (HiCanu + purge_dups + YaHS); k-mer completeness (68 %) as principal limitation incl. Reviewer B's expectation model and its caveats; Methods list every departure from the preregistration. Checked against code: adjacent-pair test 2,000 permutations; windowed scan null = 99th percentile of the genome-wide maximum over 300 within-chromosome permutations; power model 20 replicates per (δ, φ); paralog rule 5 percentage points. Facts corrected from sources: Bonn plant collected in Braunschweig (grown from a single root; authors read the 1.1 Gb as "two haplophases"); MH1 = 567 Mb (Hirabayashi 2026 Table 1); DToL sex MALE from BioSample SAMEA7521824 (ENA), not from the note (which also gives 44.92 Gb HiFi = "42-fold", i.e. a ~1.07-Gb reference size).

**Final novelty search (web_search, 2026-10-24).** Queries: "Urtica dioica autotetraploid subgenome divergence genome"; "stinging nettle tetraploid origin allopolyploid genome assembly 2026"; "Urtica dioica polyploid origin genomic evidence haplotype"; "Urtica dioica tetraploid k-mer Smudgeplot GenomeScope autopolyploid". Results: Hirabayashi 2025/2026, Christenhusz 2025, Wolff 2026, Rejlová 2019, Cronk 2016, a 2026 Urticaceae plastome paper (Frontiers), GenomeScope2 papers — none reports a genome-wide auto/allo test or a subgenome-divergence measurement for U. dioica. The Bonn preprint states the 1.1 Gb "indicates that two haplophases are represented" — i.e. its authors read the plant as diploid. No racing conflict for H1; SDR work (UBC group's stated next step) was never a target here.

**Code packaging.** code/prep_proteomes.py reconstructed from the 2026-10-23 inline cells and verified to regenerate results/synteny/{UDI,CAN,MOR,UUR,PJU}.faa and *_genes.csv byte-identically (MD5). code/copy_structure_v2.py written and verified to regenerate results/tables/copy_structure_v2.json identically (13 min). code/sensitivity_tables.py regenerates the two small tables identically. Snakefile: rules h1_v2, copy_structure, copy_structure_v2, figures added; prep_proteomes.py now exists; trim-0 run documented as optional. make_figures.py: Fig 1 uses h1v2 chromosome-level CIs, Fig 2a includes pooled baselines, Fig 5 added. data/PROVENANCE.md completed from download_log.tsv (33 files, 25 MD5-verified, 8 without server MD5 = Ensembl/assembly reports) plus supplementary and literature fetches.

**Documents.** paper/REVIEWER_RESPONSES.md (point-by-point to A1–A4, B1–B17, B-F1/F2); paper/CLAIM.md (final wording: "no divergent second subgenome within the deposited sequence"; excluded (δ, φ) region); paper/CLAIM_LEDGER.md (23 rows, all keyed to h1v2/v2 files); SUMMARY.md; VALIDATION.md (item 5 now the first-priority read-based test with expected outcomes under each hypothesis); paper/AI_DISCLOSURE.md (review process, what was not done, human checks required).

**Not done / open.** Read-based k-mer test (VALIDATION.md 5); Rejlová 2021 full text and PhD thesis unread; Reviewer specialist criteria (platform UI setting) not changed; raw FASTA/BLAST DBs (~14 GB) left on disk pending user decision.


## 2026-10-24T14:02:20+01:00 — Accuracy pass on title and claims (user request)

Read Wolff et al. 2026 full text (bioRxiv PDF via Unpaywall; 27 pp). Facts: 53.48 Gb ONT R9 + 73.17 Gb R10; selected assembly hifiasm on R10, 1.1 Gb, 491 contigs, N50 40.67 Mb; BUSCO embryophyta_odb12 C 97.5 % [S 16.7 %, D 80.8 %]; "high percentage of duplicated BUSCO genes indicates that haplophases have been separated"; coverage of all R9 reads mapped to the assembly peaks at ≈ 25× and ≈ 50×, mean ≈ 40×; MGSE genome size ≈ 1.3 Gb from 53.48 Gb; NextDenovo run with expected genome size 1.6 Gb; assembly "substantially smaller than the previous estimate"; plant XX-0-BONN-49798 from Braunschweig, grown from a single root (no pooling). Arithmetic: mean 40× × 1.1 Gb ≈ 44 Gb mapped; diploid model predicts hap-separated/collapsed peaks at 44/(2·Cx) and 44/Cx = 37×/73× for Cx = 0.6 Gb; tetraploid model predicts 18×/37×/73×. Observed 25×/50× is ~1/3 off in opposite directions under the two models → the depth figure does not discriminate as published. BUSCO D 80.8 % with two copies dominant is expected under both models. MGSE is doubled by two haplophases in the reference under both models. Conclusion: Claim 2 downgraded from "probable tetraploid" to "ploidy unresolved; tetraploidy possible". Title changed to "No second subgenome more divergent than diploid alleles in the deposited genome of a tetraploid stinging nettle (Urtica dioica): a preregistered genome-wide test against haplotype-resolved diploid assemblies" so that it states exactly the tested quantity, the sample (one deposited genome) and the comparator. "first genome-wide" → "to our knowledge the first". Rejlová 2021 full text: closed access (Unpaywall/S2/PMC/CrossRef-TDM all negative); still cited from abstract.


## 2026-10-24T15:39:47+01:00 — Rejlová 2021 full text, author block, reference style, submission PDF, deposit bundle

Rejlová et al. 2021 full text (user-supplied PDF, 16 pp): Hyb-Seq (Sondovač-designed LCN probes), ASTRAL species tree from 249 gene trees — European subsp. dioica tetraploids plus diploid subspp. kurdistanica/subinermis/pubescens/sondenii collapse to a polytomy sister to seven Asian species; 427,633 SNP sites → 26,052 biallelic → 2,135 independent; PCoA/fastStructure: no structure among diploids; Middle-East tetraploids separate. No allele phasing, no subgenome analysis (consensus per individual). Quotes used: "probably of autopolyploid origin" (Po-basin tetraploids with subsp. pubescens), "The evolutionary history of the remaining European tetraploids remains unresolved", "either the auto- or allopolyploid origin of tetraploids". Crown clade 2.4–2.0 Mya (Wu et al. 2018, cited there). Introduction and Discussion updated accordingly; our framing ("first genome-wide assembly-based test"; their design pooled alleles) confirmed.

Author block: Edward Mason, University of Oxford, self-funded, edward.mason@ell.ox.ac.uk; AI use declared in the author block and in a full AI disclosure section inside the PDF. Reference style: author–year in text (Harvard-type, as used by Botanical Journal of the Linnean Society / New Phytologist), alphabetical list; consortium authors of the three Wellcome Open Research notes restored from CrossRef (9–10 authors → "Christenhusz et al."), 2024a/2024b disambiguated. Abstract shortened 510 → 371 words. Submission PDF built with code/build_pdf.py (markdown → HTML → xhtml2pdf; DejaVu Sans embedded for Unicode; figures embedded after their legends; AI disclosure and supplementary list as sections): 17 pages, 1.98 MB, zero render errors. Deposit tree deposit/urtica-tetraploid-origin (118 files, 120 MB; largest gene_distances.csv 33 MB — under GitHub's 100 MB limit) and bioRxiv supplementary ZIP (29.5 MB).


## 2026-10-24T16:16:11+01:00 — Abstract length vs bioRxiv norms

bioRxiv publishes no abstract word limit in its FAQ (web search 2026-10-24; FAQ says submitted abstract metadata are later replaced by text extracted from the PDF). Empirical norm: api.biorxiv.org details for 2026-09-01..30, first 6,000 of 6,874 records (cap), deduplicated to 5,693 preprints with ≥ 30-word abstracts; word count = whitespace tokens containing a letter or digit (same rule applied to ours). Median 204 words (IQR 169–245, P90 278; 5.8 % > 300); genomics median 203 (n 201), evolutionary biology 216 (n 246), plant biology 214 (n 188), genetics 200 (n 129). Previous abstract 343 words = 98th–99th percentile. Rewritten to 266 words (≈ 80th–90th percentile), numbers regenerated from h1v2_result_*.json and manuscript_numbers.json; wording tightened ("not been tested genome-wide" rather than "with genome sequences", since Rejlová 2021 used Hyb-Seq; "median divergence … times"; "≥ 20 % under the preregistered allotetraploid criterion"). Previous version kept in notes/abstract_v2_371words.md. Table: results/tables/biorxiv_abstract_lengths_2026-09.csv.
