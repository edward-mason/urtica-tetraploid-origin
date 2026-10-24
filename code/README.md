# code/ — reproducible workflow

Windows-native (no WSL). Requirements: Python 3.11 environment from `environment.yml`; NCBI BLAST+ 2.17.0 x64 Windows build unpacked into `code/bin/` (download: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-win64.tar.gz; MD5 in the same directory).

Run everything: `snakemake -s code/Snakefile --cores 8 all` (≈ 30 h on 4 cores/8 threads; peak disk ≈ 14 GB; one genome is decompressed and indexed at a time and the plain FASTA is deleted after indexing).

Steps and scripts
1. `download.py handoff/download_plan.json data/genomes` — serial download with MD5 verification (NCBI FTP, Ensembl Rapid Release); `data/PROVENANCE.md` lists bonndata/figshare/Europe PMC items fetched separately.
2. `gene_divergence.py prepare` — longest RefSeq CDS per gene (MH2) → `results/divergence/query_cds.fa`, `query_index.csv`.
3. `gene_divergence.py blast <T>` for T in MH2 FH1 FH2 MH1 DP DA UU BN — makeblastdb + dc-megablast (evalue 1e-20, max_target_seqs 5, max_hsps 40).
4. `gene_divergence.py summarise <T>` — locus clustering, single-copy ratio, query-anchored sequences (trim 10).
5. `h1_distances.py <targets>` — pairwise p-distances per gene, parental labels (Hirabayashi 2026 Table 2) → `gene_distances.csv`.
6. `h1_test.py DP|DA|BN` — preregistered H1 statistics → `h1_result_*.json`, `h1_genes_*.csv`.
6b. `h1_test_v2.py DP|DA` — post-review statistics (paralog screen, chromosome-level bootstrap CIs, adjacent-pair clustering permutation, windowed paired t−b scan, Hartigan dip test, (δ, φ) power model, chr 3/5/6 symmetry check) → `h1v2_result_*.json`, `h1v2_genes_*.csv`. These are the numbers reported in the manuscript.
7. `copy_structure.py` — loci per gene at ≥97 % identity / ≥80 % coverage → `copy_structure.json`; `copy_structure_v2.py` — pooled unphased-diploid (FH1+FH2, MH1+MH2) and pooled-tetraploid (DP+DA) baselines and Bonn contig-length context → `results/tables/copy_structure_v2.json`.
8. `h3_test.py female|male` — preregistered H3 → `h3_result_*.json`, `h3_perchr_*.csv` (needs `data/annot/centromere_points*.csv`, built from the published supplementary tables as described in notes/notebook.md).
9. `prep_proteomes.py` (longest protein per gene + gene positions for UDI/UUR/PJU/CAN/MOR → `results/synteny/`), then `h2_synteny.py blast|blocks <SP>` and `h2_test.py` — RBH synteny and the preregistered H2 rule → `results/h2_result.json`.
Figures: `make_figures.py` (re-creates Fig1–Fig5 from the result files).
Sensitivity analyses: trim-0 re-run (`NETTLE_TRIM=0 NETTLE_OUTDIR=results/divergence_trim0`, then summarise/distances/h1 for DP) → `results/tables/h1_trim0_sensitivity.json`; residual t > 5 % genes → `results/tables/h1_residual_t_gt5_sensitivity.json`; second-locus identity classes → `results/tables/second_locus_identity_classes.json` (`sensitivity_tables.py`).
Independent reproduction scripts (Reviewer A): `paper/reviews/reviewerA_tier1.py`, `reviewerA_tier2.py`.
Seeds: 20261023 throughout. Preregistration: `notes/preregistration.md` (SHA-256 in `notes/preregistration.sha256`).
