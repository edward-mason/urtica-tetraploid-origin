# data_PROVENANCE_reads.md — Public sequencing runs for *Urtica* (NCBI taxid 3500 subtree)

Query date: 2026-10-22. All counts below come from API responses retrieved on this date; files in `reads_inv/` hold the raw responses.

## 1. Sources and queries

| Source | Query | Result |
|---|---|---|
| ENA portal API `search?result=read_run&query=tax_tree(3500)&limit=0` (35 fields) | all runs in taxid 3500 subtree | 329 runs, 329 experiments, 228 samples, 34 studies |
| NCBI E-utilities `esearch db=sra term=txid3500[Organism]` | cross-check | 329 SRA experiments (exact match with ENA) |
| NCBI `esearch db=sra term=txid3501[Organism]` | *U. dioica* incl. subspecies | 237 experiments; ENA: 237 runs with scientific_name starting "Urtica dioica" |
| ENA browser XML for all 34 studies and all 228 samples | publication xrefs, sample attributes | 34/34 and 228/228 retrieved |
| NCBI `esummary`/`efetch db=bioproject` for 34 projects | descriptions, publications | 34/34; 0 `<Publication>` elements |
| Europe PMC search for each study accession (primary + secondary) | full-text accession mentions | candidate papers for 23 studies (text mentions, not necessarily the primary paper) |
| NCBI Datasets v2 `genome/taxon/3500/dataset_report` | assemblies + BioSample links | 8 assemblies |
| NCBI `esearch db=gds` (Urtica[Organism], nettle) | GEO | 2 relevant GSE series (GSE185155, GSE81702) |
| NCBI nuccore (tsa filters), ENA `result=tsa_set`, `result=analysis`, EVA species list | TSA / VCF / analyses | 0 TSA records; 0 ENA analyses; *Urtica* absent from EVA |
| Dryad API, Zenodo API (`q="Urtica dioica"`, first 100 of 1,872 text hits) | processed data | no sequencing-derived *Urtica* datasets found (Zenodo scan partial) |

Notes on the ENA record: 329 runs map 1:1 to experiments. 25 runs have `read_count`=0 and `base_count`=0 in ENA metadata (24 ONT native-format WGA runs in PRJEB63450 and Hi-C run SRR31536218), so base totals for those are under-reported. The *Urtica dioica* subspecies in the ENA set: subsp. dioica (54 runs), subinermis (18), pubescens (9), sondenii (4), kurdistanica (2), var. angustifolia (3; NCBI taxid 1399755, an East Asian taxon treated by some authors as the separate species *U. angustifolia* — these three runs are counted in the 237 by their ENA scientific_name but are not Eurasian *U. dioica* s.str.), plus 147 runs labelled only "Urtica dioica".

## 2. Per-category totals

Classification rule: `RNA-Seq`→RNA-seq; `Hi-C`→Hi-C; `Bisulfite-Seq`→methylation; `WGS`/`WGA` on PACBIO_SMRT or OXFORD_NANOPORE→WGS long-read; `WGS` on ILLUMINA/DNBSEQ→WGS short-read (includes genome skims and plastome projects); `AMPLICON`, `Targeted-Capture`, `OTHER` (= HybSeq per experiment titles)→amplicon/other.

### 2a. Whole genus (taxid 3500 subtree)
| category | runs | samples | studies | Gb | runs_with_zero_or_missing_bases |
|---|---|---|---|---|---|
| amplicon/other (target capture, HybSeq, amplicon) | 125 | 123 | 10 | 29.1 | 0 |
| RNA-seq | 118 | 58 | 9 | 295.39 | 0 |
| WGS short-read | 45 | 39 | 14 | 313.13 | 0 |
| WGS long-read | 31 | 6 | 5 | 212.91 | 24 |
| Hi-C | 5 | 4 | 4 | 439.58 | 1 |
| methylation (bisulfite/EM-seq) | 5 | 5 | 1 | 0.57 | 0 |

### 2b. *Urtica dioica* s.l. only (scientific_name starts with "Urtica dioica")
| category | runs | samples | studies | Gb | runs_with_zero_or_missing_bases |
|---|---|---|---|---|---|
| RNA-seq | 98 | 38 | 6 | 117.68 | 0 |
| amplicon/other (target capture, HybSeq, amplicon) | 85 | 84 | 4 | 23.33 | 0 |
| WGS long-read | 30 | 5 | 4 | 194.05 | 24 |
| WGS short-read | 15 | 13 | 7 | 159.19 | 0 |
| methylation (bisulfite/EM-seq) | 5 | 5 | 1 | 0.57 | 0 |
| Hi-C | 4 | 3 | 3 | 316.03 | 1 |

Long-read detail (*U. dioica*): PacBio HiFi Revio — SRR31536217 (clone 11-4, 66.3 Gb) and SRR35053641 (clone 48M, 77.5 Gb); PacBio Sequel II/IIe — ERR10662027, ERR10662028 (drUrtDioi1, 43.4 Gb total); ONT MinION — SRR31035787 (clone 11-4, 5.3 Gb) plus 24 whole-genome-amplified MinION runs in PRJEB63450 (sizes not reported in ENA). No methylation data other than the 5 bisulfite runs of PRJNA322330 (0.57 Gb total; leaf; no origin metadata). No EM-seq, no ONT methylation calls deposited.

## 3. Studies (34; sorted by total bases) — full table in `studies.csv`
| study_accession | n_runs | n_samples | total_Gb | categories | species | first_public | publication |
|---|---|---|---|---|---|---|---|
| PRJEB58082 | 7 | 2 | 345.76 | Hi-C; RNA-seq; WGS long-read | Urtica dioica | 2022-12-13 | VERIFIED 10.12688/wellcomeopenres.23399.1 |
| PRJEB47304 | 12 | 7 | 270.8 | Hi-C; RNA-seq; WGS long-read; WGS short-read | Urtica urens | 2021-09-10 | none found |
| PRJNA663211 | 5 | 1 | 204.97 | Hi-C; WGS long-read; WGS short-read | Urtica dioica | 2023-10-02 | VERIFIED 10.3390/plants14010124 |
| PRJNA1308592 | 2 | 1 | 109.88 | Hi-C; WGS long-read | Urtica dioica subsp. dioica | 2025-08-24 | VERIFIED 10.1038/s41597-026-06960-7 |
| PRJNA1139937 | 12 | 12 | 91.62 | RNA-seq | Urtica cannabina | 2024-07-29 | (mention only) 10.3389/fpls.2025.1598628 (2025) The short-term effect  |
| PRJDB41448 | 22 | 22 | 69.33 | WGS short-read | Urtica ardens; Urtica atrichocaulis; Urtica chamaedryoides; Urtica dio | 2026-06-02 | none found |
| PRJEB63450 | 29 | 7 | 40.15 | RNA-seq; WGS long-read | Urtica dioica | 2024-12-31 | none found |
| PRJNA592832 | 2 | 2 | 38.15 | RNA-seq | Urtica ferox; Urtica incisa | 2021-01-04 | VERIFIED 10.1126/sciadv.abb8828 |
| PRJNA518071 | 76 | 16 | 23.77 | RNA-seq | Urtica dioica | 2020-02-26 | VERIFIED 10.1002/pld3.151 |
| PRJNA602985 | 87 | 87 | 23.71 | amplicon/other (target capture, HybSeq, amplicon) | Urtica atrovirens; Urtica cypria; Urtica dioica subsp. dioica; Urtica  | 2020-02-29 | (mention only) 10.3390/ijms21113853 (2020) Identification of Callose S |
| PRJNA1237398 | 2 | 2 | 20.76 | RNA-seq | Urtica dioica | 2026-04-19 | (mention only) 10.1099/jgv.0.002273 (2026) Variability and population  |
| PRJNA767836 | 12 | 12 | 13.07 | RNA-seq | Urtica dioica | 2021-11-17 | VERIFIED 10.3390/ijms222212319 |
| PRJEB55612 | 1 | 1 | 7.59 | WGS short-read | Urtica membranacea | 2022-09-06 | (mention only) 10.1111/tpj.70459 (2025) Progestogens and androgens inf |
| PRJEB43865 | 5 | 5 | 4.72 | WGS short-read | Urtica dioica; Urtica dioica subsp. dioica; Urtica dioica subsp. sonde | 2021-08-17 | (mention only) 10.3390/ijms26189240 (2025) Integration of Repeatome an |
| PRJEB82787 | 3 | 2 | 4.6 | WGS short-read | Urtica dioica; Urtica urens | 2024-11-27 | (mention only) 10.3389/fpls.2026.1821377 (2026) Comparative chloroplas |
| PRJNA749813 | 1 | 1 | 4.51 | WGS short-read | Urtica dioica var. angustifolia | 2021-07-28 | (mention only) 10.1080/23802359.2022.2057246 (2023) The complete chlor |
| PRJNA782209 | 1 | 1 | 4.42 | WGS short-read | Urtica fissa | 2021-11-22 | (mention only) 10.1080/23802359.2022.2080017 (2022) The complete chlor |
| PRJEB21674 | 1 | 1 | 2.75 | RNA-seq | Urtica dioica | 2017-07-23 | (mention only) 10.3390/v15040840 (2023) Expanding the Repertoire of th |
| PRJDB33230 | 1 | 1 | 1.61 | WGS short-read | Urtica radicans | 2025-09-28 | none found |
| PRJEB88382 | 4 | 2 | 1.45 | amplicon/other (target capture, HybSeq, amplicon) | Urtica dioica; Urtica urens | 2026-08-25 | (mention only) 10.1111/1755-0998.70202 (2026) Simultaneous Inference o |
| PRJNA1022032 | 11 | 11 | 1.25 | amplicon/other (target capture, HybSeq, amplicon) | Urtica gracilenta; Urtica hyperborea; Urtica kioviensis; Urtica laetev | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJEB30946 | 1 | 1 | 1.02 | WGS short-read | Urtica dioica | 2019-02-08 | none found |
| PRJNA1022023 | 8 | 8 | 0.9 | amplicon/other (target capture, HybSeq, amplicon) | Urtica atrovirens; Urtica cannabina; Urtica chamaedryoides; Urtica cir | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJEB50550 | 1 | 1 | 0.72 | WGS short-read | Urtica gracilis subsp. holosericea | 2022-09-20 | none found |
| PRJNA322330 | 5 | 5 | 0.57 | methylation (bisulfite/EM-seq) | Urtica dioica | 2017-03-17 | VERIFIED 10.1186/s13059-017-1195-1 |
| PRJNA325670 | 1 | 1 | 0.51 | WGS short-read | Urtica dioica | 2023-04-02 | VERIFIED 10.1055/s-0043-113449 |
| PRJNA1022141 | 6 | 6 | 0.44 | amplicon/other (target capture, HybSeq, amplicon) | Urtica echinata; Urtica glomeruliflora; Urtica longispica; Urtica macb | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJNA1022138 | 3 | 3 | 0.43 | amplicon/other (target capture, HybSeq, amplicon) | Urtica atrichocaulis; Urtica mairei; Urtica praetermissa | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJNA1022323 | 2 | 2 | 0.43 | amplicon/other (target capture, HybSeq, amplicon) | Urtica ardens; Urtica massaica | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJNA1022147 | 1 | 1 | 0.28 | amplicon/other (target capture, HybSeq, amplicon) | Urtica gracilis subsp. mollis | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJEB37665 | 1 | 1 | 0.22 | WGS short-read | Urtica urens | 2020-07-13 | none found |
| PRJNA1022025 | 2 | 2 | 0.2 | amplicon/other (target capture, HybSeq, amplicon) | Urtica andicola; Urtica gracilis subsp. aquatica | 2023-11-10 | (mention only) 10.1093/jhered/esag028 (2026) The origin of the octoplo |
| PRJEB114749 | 1 | 1 | 0.08 | WGS short-read | Urtica urens | 2026-06-19 | none found |
| PRJNA649768 | 1 | 1 | 0.0 | amplicon/other (target capture, HybSeq, amplicon) | Urtica dioica | 2020-08-01 | (mention only) 10.1002/aps3.11405 (2021) Examining the utility of DNA  |

Publication linkage: only PRJNA322330 (PMID 28457232, doi:10.1186/s13059-017-1195-1), PRJNA325670 (PMID 28651291, doi:10.1055/s-0043-113449) and PRJNA767836 (PMID 34830202, doi:10.3390/ijms222212319) carry explicit PubMed xrefs in ENA study XML. For the other studies `studies.csv` lists Europe PMC full-text mentions; the ones I consider firmly tied by content (data-availability statements read in full text) are: PRJNA663211 → doi:10.3390/plants14010124; PRJNA1308592 → doi:10.1038/s41597-026-06960-7; PRJNA518071 → doi:10.1002/pld3.151; PRJNA592832 → doi:10.1126/sciadv.abb8828; PRJEB58082 → DToL genome note doi:10.12688/wellcomeopenres.23399.1 (Crossref abstract; the note does not appear in Europe PMC accession search). PRJNA602985 (HybSeq, Institute of Botany CAS) has no confirmed publication; the Europe PMC hit (doi:10.3390/ijms21113853) is a text mention only.

## 4. Cytotype per BioSample — full table in `biosamples_cytotype.csv`

No BioSample in the set carries a `ploidy`, chromosome-number or karyotype attribute (regex over all attributes: 1 hit, which was a false positive on "sub_species"). Verdicts therefore rest on linked assemblies and publications.

Totals (228 BioSamples): diploid = 2, tetraploid = 2, unknown = 136 (85 with a weak subspecies hint, 51 with none), n/a non-dioica species = 88. **All 136 "unknown" samples are UNUSABLE for ploidy-sensitive analyses.**

### 4a. Samples with a verdict (evidence a = explicit statement; b = same individual as a karyotyped/chromosome-level assembly)
| sample_accession | study_accession | sci_name | sex | geo | strategies | cytotype_verdict | evidence_level | evidence |
|---|---|---|---|---|---|---|---|---|
| SAMEA7521894 | PRJEB58082 | Urtica dioica | MALE | United Kingdom | RNA-Seq/WGS | tetraploid | a+b | Tissue sample "derived from" SAMEA7521824 (ToLID drUrtDioi1). DToL genome note (doi:10.12688/wellcomeopenres.23399.1) abstract: 618.00 Mb assembly, 13 pseudomolecules "representing one subgenome of the tetraploid genome". Corroborated by Sci Data 2026: "A haplotype from a tetraploid individual ... was also subsequently generated". Assembly GCA_964188135.1 (+ alt GCA_964188125.1). |
| SAMEA7521898 | PRJEB58082 | Urtica dioica | MALE | United Kingdom | Hi-C/RNA-Seq | tetraploid | a+b | Tissue sample "derived from" SAMEA7521824 (ToLID drUrtDioi1). DToL genome note (doi:10.12688/wellcomeopenres.23399.1) abstract: 618.00 Mb assembly, 13 pseudomolecules "representing one subgenome of the tetraploid genome". Corroborated by Sci Data 2026: "A haplotype from a tetraploid individual ... was also subsequently generated". Assembly GCA_964188135.1 (+ alt GCA_964188125.1). |
| SAMN16122873 | PRJNA663211 | Urtica dioica | female | Romania: River Jiu, north of Rovinari. | Hi-C/WGS | diploid | a+b | Plants 2025 (doi:10.3390/plants14010124): "chromosome-level, phased genome assembly for a diploid female clone of Urtica dioica from Romania"; flow cytometry 2C=1.26 pg -> 1C=616 Mb; 2n=26. Same BioSample as assemblies GCA_048565285.1/GCA_048565295.1 (13 chromosomes, 575/521 Mb). |
| SAMN50702769 | PRJNA1308592 | Urtica dioica subsp. dioica | male | Canada: Vancouver | Hi-C/WGS | diploid | a+b | Sci Data 2026 (doi:10.1038/s41597-026-06960-7): "haplotype-resolved sequence of a diploid male plant" (clone 48M, offspring of diploid 11-4 x diploid 7-5). Same BioSample as GCA_054948195.1/GCA_054948205.1 (13 chromosomes, 567/534 Mb; RefSeq GCF_054948205.1). |

Interpretation: the two Cronk-lab individuals (clone 11-4, Romania; clone 48M = 11-4 × Bulgarian male 7-5) are diploid 2n=26, 1C≈616 Mb by flow cytometry, with 13-chromosome assemblies of 521–575 Mb per haplotype. The Darwin Tree of Life specimen drUrtDioi1 (Kew, Surrey, UK; male) is tetraploid: its 618 Mb, 13-pseudomolecule assembly GCA_964188135.1 is explicitly "one subgenome of the tetraploid genome". Consequently the DToL HiFi/Hi-C/RNA-seq reads (PRJEB58082) come from a tetraploid individual while the assembly covers one subgenome — any analysis contrasting the two cytotypes must treat PRJEB58082 reads and the drUrtDioi1 assembly differently. Caveat: the Cronk-lab diploids are labelled "subsp. dioica" in INSDC, which contradicts the taxonomic convention that subsp. dioica is the tetraploid; subspecies labels are therefore only weak hints.

### 4b. *U. dioica* samples by study and verdict
| study_accession | cytotype_verdict | evidence_level | n_samples | weak_hint |
|---|---|---|---|---|
| PRJDB41448 | unknown | c (weak) | 2 | labelled subsp. dioica -> nominate subspecies is the widespr |
| PRJDB41448 | unknown | none | 1 |  |
| PRJEB21674 | unknown | none | 1 |  |
| PRJEB30946 | unknown | none | 1 |  |
| PRJEB43865 | unknown | c (weak) | 3 | labelled subsp. dioica -> nominate subspecies is the widespr; labelled subsp. sondenii -> reported diploid (2n=26) subsp |
| PRJEB43865 | unknown | none | 1 |  |
| PRJEB58082 | tetraploid | a+b | 2 |  |
| PRJEB63450 | unknown | none | 7 |  |
| PRJEB82787 | unknown | none | 1 |  |
| PRJEB88382 | unknown | none | 1 |  |
| PRJNA1022023 | unknown | none | 1 |  |
| PRJNA1237398 | unknown | none | 2 |  |
| PRJNA1308592 | diploid | a+b | 1 |  |
| PRJNA322330 | unknown | none | 5 |  |
| PRJNA325670 | unknown | none | 1 |  |
| PRJNA518071 | unknown | none | 16 |  |
| PRJNA602985 | unknown | c (weak) | 80 | labelled subsp. dioica -> nominate subspecies is the widespr; labelled subsp. kurdistanica -> reported diploid (2n=26) s |
| PRJNA649768 | unknown | none | 1 |  |
| PRJNA663211 | diploid | a+b | 1 |  |
| PRJNA767836 | unknown | none | 12 |  |

Weak-hint basis (evidence c): Rejlová et al. 2019, PLoS One 14:e0218389 (doi:10.1371/journal.pone.0218389) — widespread tetraploid cytotype 2n=4x=52 = subsp. dioica; diploid subspecies 2n=2x=26 = subsp. subinermis, pubescens, sondenii, kurdistanica, cypria (subinermis also has a reported tetraploid cytotype); triploids/pentaploids occur sporadically. Cronk et al. 2016 (doi:10.3897/BDJ.4.e10003): 5/29 transect plants diploid (Bulgaria, Romania, S. Europe), 2C diploid 1.31–1.35 pg, tetraploid 2.36–2.59 pg. Neither paper was used to assign a firm verdict to any BioSample.

## 5. Published processed outputs usable without read alignment
| Source | Content | File(s) | Size | URL | Verification status |
|---|---|---|---|---|---|
| GEO GSE185155 (PRJNA767836, callogenesis RNA-seq, clone 13) | De novo transcriptome assembly | GSE185155_Callinettle_assembly.fa.gz | 9.5 MB | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185155/suppl/ | Listed via FTP; not downloaded |
| GEO GSE185155 | Per-sample "Total gene reads" + RPKM for 12 samples (0D/1D/10D/34D x3), 3,143 contigs only (CLC export; appears DE-filtered, not a full matrix) + ANOVA stats + annotation | GSE185155_Processed_RNA-Seq_data_nettle_callogenesis.txt.gz | 1.1 MB (3,143 rows, 105 cols) | same dir | Downloaded and header verified |
| GEO GSE185155 | ANOVA statistics sub-experiment (xlsx); blast2go annotation | GSE185155_0D_vs._1D_vs._10D_vs._34D_ANOVA_stat_subexperiment_GEO.xlsx; GSE185155_blast2go_export_...txt.gz | 29 MB; 3.4 MB | same dir | Listed via FTP |
| GEO GSE81702 (PRJNA322330, Bewick et al. 2017 Genome Biol, doi:10.1186/s13059-017-1195-1) | Per-cytosine methylation calls (chromosome, position, strand, context, methylated reads, total reads, binomial test) for 5 U. dioica leaf replicates GSM2171637-GSM2171641, inside a multi-species tarball (150 samples) | GSE81702_processed_data_files.tar.gz | 38 MB (exceeds 20 MB limit; not downloaded) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81702/suppl/ | Metafile downloaded (column definitions). Reference coordinate system for U. dioica NOT stated in metadata (paper says browsers at Plant Methylation DB, schmitzlab; not reachable) |
| NCBI RefSeq annotation GCF_054948205.1 RS_2026_09 (clone 48M H2, diploid male) | Gene models: GFF/GTF, 30,210 XM_ mRNAs (nuccore count), proteins, feature table; plus list of 70 SRA RNA-seq runs NCBI aligned, with per-run alignment rates | *_genomic.gff.gz 7.3 MB; *_protein.faa.gz 6.5 MB; *_rna.fna.gz 14 MB; *_feature_table.txt.gz 2.3 MB; RS_2026_09_rnaseq_runs.txt 14 KB; rnaseq_alignment_summary.txt 4.4 KB | see sizes | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/ | Two small txt files downloaded; files dated 2026-10-20 |
| Figshare project 230981 (Cronk lab; cited in doi:10.3390/plants14010124 and doi:10.1038/s41597-026-06960-7) | Gene + TE annotation GFF for clone 11-4 H1/H2 and clone 48M MH1/MH2, assemblies, supplementary tables | not enumerated | unknown | https://figshare.com/projects/A_high-quality_phased_genome_assembly_of_stinging_nettle_Urtica_dioica_ssp_dioica/230981 | api.figshare.com blocked by allowlist; figshare.com page is JS-only (HTTP 202, empty). File list NOT verified |
| GenBank assemblies (NCBI Datasets, taxon 3500) | 8 assemblies: 3 U. dioica individuals x haplotypes (GCA_964188135.1/125.1 drUrtDioi1; GCA_048565285.1/295.1 clone 11-4; GCA_054948195.1/205.1 clone 48M) + U. urens drUrtUren1 (GCA_958296335.1/315.1) | see reads_inv/urtica_assemblies_ncbi.csv | n/a | https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=3500 | Metadata verified |
| Frontiers Plant Sci 2025 doi:10.3389/fpls.2025.1598628 (PRJNA1139937, U. cannabina root RNA-seq) | Supplementary Table 7 "The table of gene reads count" | Supplementary Table 7 | unknown | publisher supplementary (frontiersin.org not reachable) | Not verified |
| 1KP sample WKCY (PRJEB21674, ERR2040431, U. dioica leaf) | Assembled transcriptome hosted by the 1KP project (outside INSDC; no TSA record exists) | n/a | unknown | 1KP data portal / CNGBdb | Not reachable from sandbox; NOT verified |
| bioRxiv 2026 preprint doi:10.64898/2026.05.15.725508 (Univ. Bonn) | 1.1 Gbp ONT-based U. dioica assembly, N50 40.7 Mb, with annotation (abstract only) | n/a | n/a | bioRxiv | Data location not stated in abstract; probable but UNVERIFIED link to PRJEB63450 (Germany, ONT + 5 tissue RNA-seq); no INSDC assembly found for it as of query date |

Negative results: no VCF for any *Urtica* in ENA (`result=analysis` returns 0 rows) or EVA (genus absent from species list); no TSA (GxxX) records for *Urtica* in NCBI nuccore or ENA `tsa_set`; GEO has no *Urtica* count matrix other than the partial GSE185155 table; no methylation calls other than GSE81702.

RefSeq RS_2026_09 RNA-seq alignment rates (NCBI's own alignment to the diploid clone-48M genome; per-run values are in `runs.csv`, column `refseq_RS_2026_09_alignment_rate_pct`): PRJEB21674 1KP leaf 90 %; PRJNA767836 callus 74–81 % (n=12); PRJNA518071 stem tissues 63–82 % (n=50); PRJEB58082 DToL tetraploid leaf 71–83 % (n=2); PRJEB63450 German tissue panel 28–86 % (n=5; leaf 38, stem 28, seed 38, seedling 72, root 86). These are the only alignment-derived numbers available without local alignment.

## 6. Geography, herbarium and sex-labelled material

Continent of collection (BioSample `geo_loc_name`; 54 samples have none): Europe 145, Asia 25, North America 4. The HybSeq panel PRJNA602985 (87 samples, 1 run each, leaf, mature plants) spans 25 countries: Italy 14, Czech Republic 14, Turkey 8, Slovakia 5, Sweden 5, Germany 5, Finland 4, France 4, Austria 3, Iran 3, Norway 3, Romania 2, Slovenia 2, Greece 2, Poland 2, Estonia 2, Bosnia and Herzegovina 1, Spain 1, Cyprus 1, Croatia 1, Latvia 1, Lithuania 1, Bulgaria 1, Georgia 1, Canada 1 — the only population-scale geographic sampling of the *U. dioica* complex, but it is target-capture data (not usable for genome-wide ploidy or expression work). Asian *U. dioica* s.l. material: Turkey (8, incl. subsp. kurdistanica ×2), Iran (3), Georgia (1), all in PRJNA602985. North American *U. dioica*: SAMN27355170 (USA: Maryland, plastome skim), SAMN15679854 (USA: SD herbarium, amplicon barcode), SAMN50702769 (grown in Vancouver but of Romanian × Bulgarian parentage), SAMN13915224 (*U. gracilis*, Ontario).

Herbarium samples (9):
| sample_accession | study_accession | sci_name | geo | strategies |
|---|---|---|---|---|
| SAMEA117369864 | PRJEB82787 | Urtica urens | France | WGS |
| SAMEA117370657 | PRJEB82787 | Urtica dioica | Switzerland | WGS |
| SAMEA8200944 | PRJEB43865 | Urtica dioica subsp. sondenii | Norway | WGS |
| SAMEA8200970 | PRJEB43865 | Urtica dioica subsp. dioica | Svalbard | WGS |
| SAMEA8201449 | PRJEB43865 | Urtica dioica | Denmark | WGS |
| SAMEA8202209 | PRJEB43865 | Urtica dioica subsp. dioica | Norway | WGS |
| SAMEA8202360 | PRJEB43865 | Urtica urens | Norway | WGS |
| SAMEA9098815 | PRJEB50550 | Urtica gracilis subsp. holosericea | Norway | WGS |
| SAMN15679854 | PRJNA649768 | Urtica dioica | USA:SD,0 | AMPLICON |

Sex-labelled samples (11):
| sample_accession | study_accession | sci_name | sex | tissue | strategies | cytotype_verdict |
|---|---|---|---|---|---|---|
| SAMEA114389642 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | ROOT | RNA-Seq | n/a (non-dioica Urtica species) |
| SAMEA114389643 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | LEAF | RNA-Seq | n/a (non-dioica Urtica species) |
| SAMEA114389644 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | **OTHER PLANT TISSUE** | RNA-Seq | n/a (non-dioica Urtica species) |
| SAMEA114389645 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | STEM | RNA-Seq | n/a (non-dioica Urtica species) |
| SAMEA7521886 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | LEAF | RNA-Seq | n/a (non-dioica Urtica species) |
| SAMEA7521888 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | LEAF | Hi-C | n/a (non-dioica Urtica species) |
| SAMEA7521889 | PRJEB47304 | Urtica urens | HERMAPHRODITE MONOECIOUS | LEAF | RNA-Seq/WGS | n/a (non-dioica Urtica species) |
| SAMEA7521894 | PRJEB58082 | Urtica dioica | MALE | LEAF | RNA-Seq/WGS | tetraploid |
| SAMEA7521898 | PRJEB58082 | Urtica dioica | MALE | LEAF | Hi-C/RNA-Seq | tetraploid |
| SAMN16122873 | PRJNA663211 | Urtica dioica | female | vegetative shoots | Hi-C/WGS | diploid |
| SAMN50702769 | PRJNA1308592 | Urtica dioica subsp. dioica | male | vegetative shoots | Hi-C/WGS | diploid |

There is **no sex-contrasted expression dataset**: the only sex-labelled *U. dioica* RNA-seq runs are the two DToL male leaf libraries (ERR10659255, ERR10659257; tetraploid). The female (clone 11-4) and male (clone 48M) diploid BioSamples have only DNA libraries (HiFi, Hi-C, Illumina WGS, ONT). The 7 sex-labelled *U. urens* samples are hermaphrodite/monoecious (DToL). PRJNA518071 (76 runs; clone 13 fibre nettle; 5-week plants; tissues TOP/MID/CBOT/FBOT internode segments, fibres, callus) and PRJEB63450 (leaf, stem, seed, seedling, root) carry no sex attribute.

## 7. Gaps and limitations

1. Cytotype is established for only 4 of 140 *U. dioica* BioSamples; the two largest RNA-seq resources (PRJNA518071, PRJNA767836; both clone 13, Germany/Luxembourg) and the bisulfite set (PRJNA322330) have no ploidy evidence (the Plant Direct 2019 full text contains no ploidy/genome-size statement).
2. Twenty-four ONT WGA runs (PRJEB63450) and one Hi-C run (SRR31536218) have no read/base counts in ENA (native-format submissions); base totals for WGS long-read and Hi-C categories are under-reported.
3. api.figshare.com is outside the allowlist and figshare.com serves a JS-only page: the Figshare project 230981 file list (gene/TE GFFs for the diploid genomes) is cited from the papers but not enumerated or size-checked.
4. The 2026 bioRxiv 1.1 Gb ONT assembly (Univ. Bonn) has no INSDC assembly record yet; its association with PRJEB63450 is inferred from geography/platform and is unverified (bioRxiv full text not reachable).
5. 1KP WKCY assembled transcriptome and the Bewick 2017 Plant Methylation DB browser are hosted outside reachable domains; existence recorded, content unverified.
6. Zenodo scan covered only the first 100 of 1,872 text hits for "Urtica dioica"; Dryad (27 hits) contained no *Urtica* sequencing datasets.
7. Publication links other than the three explicit PubMed xrefs and five content-verified papers are Europe PMC text mentions and may be re-users rather than data producers.
8. PRJNA1237398 RNA-seq (Poland) targets viral RNA in nettle leaves (library_source VIRAL RNA); host reads are incidental.

## 8. Files

- `runs.csv` — 329 rows, one per run; ENA fields + category, long-read subtype, BioSample sex/tissue/geo, cytotype verdict, RefSeq alignment rate, processed-output pointer.
- `studies.csv` — 34 rows; counts, strategies, platforms, total bases, species, countries, explicit publication xrefs, content-verified publication DOIs, Europe PMC mentions, processed outputs.
- `biosamples_cytotype.csv` — 228 rows; attributes, cytotype verdict, evidence level/text, weak hint, usability flag.
- `reads_inv/` — raw ENA TSV, study/sample JSON, Europe PMC JSON, GEO table, assembly table, processed-output table, RefSeq RNA-seq run list and alignment summary, GSE185155 processed table.
