# Data provenance

Every input: resource, accession, version, URL, date, size, checksum. Companion inventories: data/PROVENANCE_assemblies.md (all public *Urtica*/Urticaceae assemblies and annotations surveyed, 2026-10-22), data/PROVENANCE_reads.md (329 SRA/ENA runs surveyed; none downloaded), data/assemblies.csv, data/annotations.csv.

## Genome assemblies and annotations (NCBI genomes FTP and Ensembl Rapid Release; downloaded 2026-10-23 by code/download.py; MD5 from the server-side md5checksums.txt / CHECKSUMS files)

Codes: FH1/FH2 = GCA_048565285.1 / GCA_048565295.1 (female 11-4, H1/H2); MH1/MH2 = GCA_054948195.1 / GCA_054948205.1 (male 48M, H1/H2; RefSeq GCF_054948205.1 annotation RS_2026_09 released 2026-10-20); DP/DA = GCA_964188135.1 / GCA_964188125.1 (DToL drUrtDioi1.1 primary/alternate; BioSample SAMEA7521824, sex MALE, RBG Kew, collected 2020-08-05); UU = GCA_958296335.1 (*U. urens* drUrtUren1.1) + Ensembl Rapid Release geneset 2024_02; PJ = GCA_954871525.1 (*Parietaria judaica* drParJuda1.1) + Ensembl 2024_02; CAN = GCF_900626175.1 (*Cannabis sativa* cs10, RefSeq); MOR = GCF_000414095.1 (*Morus notabilis*, scaffold-level; downloaded, not used).

| downloaded (local time) | code | file | URL | bytes | MD5 | verified |
|---|---|---|---|---|---|---|
| 2026-10-23T08:42:44 | FH1 | GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/565/285/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_assembly_report.txt | 157,289 | bb5c951ed4038af32b37f0d997b03c17 | True |
| 2026-10-23T08:43:32 | FH1 | GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/565/285/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_genomic.fna.gz | 161,216,736 | 9ef0011b91a0dcfbdcfa09182eaa4855 | True |
| 2026-10-23T08:43:33 | FH2 | GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/565/295/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_assembly_report.txt | 27,680 | 98ae5f2e03e7a85b42e68fad6f57ed56 | True |
| 2026-10-23T08:44:12 | FH2 | GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/565/295/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_genomic.fna.gz | 146,416,765 | 9462f874c200cb73f05c907c2cba5f48 | True |
| 2026-10-23T08:44:13 | MH1 | GCA_054948195.1_U.dioica_male_clone48M_v1_H1_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/948/195/GCA_054948195.1_U.dioica_male_clone48M_v1_H1/GCA_054948195.1_U.dioica_male_clone48M_v1_H1_assembly_report.txt | 75,179 | 8c7db3424d5b0209bdd458d6b506cd2e | True |
| 2026-10-23T08:44:51 | MH1 | GCA_054948195.1_U.dioica_male_clone48M_v1_H1_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/948/195/GCA_054948195.1_U.dioica_male_clone48M_v1_H1/GCA_054948195.1_U.dioica_male_clone48M_v1_H1_genomic.fna.gz | 159,698,164 | cf16dda2b51eddc1c487cae8f4c60436 | True |
| 2026-10-23T08:44:53 | MH2 | GCF_054948205.1_U.dioica_male_clone48M_v1_H2_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_assembly_report.txt | 9,066 | 1516120ecf2bc8939affdb670739ff64 | True |
| 2026-10-23T08:44:58 | MH2 | GCF_054948205.1_U.dioica_male_clone48M_v1_H2_cds_from_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_cds_from_genomic.fna.gz | 12,037,643 | d4918d3fd16ef67169884b3a2e01507d | True |
| 2026-10-23T08:45:25 | MH2 | GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.fna.gz | 148,816,171 | 0a9da87a99bc64f5d5c0036c886c534a | True |
| 2026-10-23T08:45:26 | MH2 | GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.gff.gz | 7,621,510 | 08130f5ac8eeda8ac5ecdc6ec2f9cbdf | True |
| 2026-10-23T08:45:28 | MH2 | GCF_054948205.1_U.dioica_male_clone48M_v1_H2_protein.faa.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/054/948/205/GCF_054948205.1_U.dioica_male_clone48M_v1_H2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_protein.faa.gz | 6,775,572 | 9007fe5d3ba3e33847e49721aa8416dd | True |
| 2026-10-23T08:45:29 | DP | GCA_964188135.1_drUrtDioi1.1_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/135/GCA_964188135.1_drUrtDioi1.1/GCA_964188135.1_drUrtDioi1.1_assembly_report.txt | 293,893 | 72a0d53c7e01b0cea23c669a48b9691b | True |
| 2026-10-23T08:46:15 | DP | GCA_964188135.1_drUrtDioi1.1_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/135/GCA_964188135.1_drUrtDioi1.1/GCA_964188135.1_drUrtDioi1.1_genomic.fna.gz | 171,876,972 | 57e9c8c940045cef43e22a6e3d17ba30 | True |
| 2026-10-23T08:46:16 | DA | GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/125/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_assembly_report.txt | 69,894 | 1a84adaadbb1fea32d882ebcf30fc2fb | True |
| 2026-10-23T08:46:50 | DA | GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/188/125/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_genomic.fna.gz | 137,527,687 | 62d75f2091c2a881fdbd5d8d2eb425a4 | True |
| 2026-10-23T08:46:51 | UU | GCA_958296335.1_drUrtUren1.1_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/958/296/335/GCA_958296335.1_drUrtUren1.1/GCA_958296335.1_drUrtUren1.1_assembly_report.txt | 9,379 | 8840842d320c6029ee7d4454b27e2a4f | True |
| 2026-10-23T08:47:18 | UU | GCA_958296335.1_drUrtUren1.1_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/958/296/335/GCA_958296335.1_drUrtUren1.1/GCA_958296335.1_drUrtUren1.1_genomic.fna.gz | 100,339,755 | 9a418c8e4c657b51cda32b75b5b0690d | True |
| 2026-10-23T08:47:19 | PJ | GCA_954871525.1_drParJuda1.1_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/954/871/525/GCA_954871525.1_drParJuda1.1/GCA_954871525.1_drParJuda1.1_assembly_report.txt | 9,478 | 4416d7b8df50f4b08b341c2a1e5007b6 | True |
| 2026-10-23T08:47:59 | PJ | GCA_954871525.1_drParJuda1.1_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/954/871/525/GCA_954871525.1_drParJuda1.1/GCA_954871525.1_drParJuda1.1_genomic.fna.gz | 157,547,546 | 2763ce67c2e97322c191bdc7ccd9cf9c | True |
| 2026-10-23T08:48:00 | MOR | GCF_000414095.1_ASM41409v2_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/414/095/GCF_000414095.1_ASM41409v2/GCF_000414095.1_ASM41409v2_assembly_report.txt | 2,769,646 | 409f86875bf77367b9acee78135f9c4b | True |
| 2026-10-23T08:48:03 | MOR | GCF_000414095.1_ASM41409v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/414/095/GCF_000414095.1_ASM41409v2/GCF_000414095.1_ASM41409v2_genomic.gff.gz | 7,406,034 | f4bd2fa58aef563fa1a94b903ee04aa0 | True |
| 2026-10-23T08:48:06 | MOR | GCF_000414095.1_ASM41409v2_protein.faa.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/414/095/GCF_000414095.1_ASM41409v2/GCF_000414095.1_ASM41409v2_protein.faa.gz | 6,891,199 | 404860d37695a5e16d8a68095f57f942 | True |
| 2026-10-23T08:48:07 | CAN | GCF_900626175.1_cs10_assembly_report.txt | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/626/175/GCF_900626175.1_cs10/GCF_900626175.1_cs10_assembly_report.txt | 21,247 | 581ee64cbca59520e5556715010e649c | True |
| 2026-10-23T08:48:14 | CAN | GCF_900626175.1_cs10_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/626/175/GCF_900626175.1_cs10/GCF_900626175.1_cs10_genomic.gff.gz | 8,006,936 | f64af0fa30c9fb8990d03b52291a9c64 | True |
| 2026-10-23T08:48:17 | CAN | GCF_900626175.1_cs10_protein.faa.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/626/175/GCF_900626175.1_cs10/GCF_900626175.1_cs10_protein.faa.gz | 7,367,821 | 69e662ed9a8896e166aa7518aea1f638 | True |
| 2026-10-23T08:48:21 | UU_ens | Urtica_urens-GCA_958296335.1-2024_02-cds.fa.gz | https://ftp.ensembl.org/pub/rapid-release/species/Urtica_urens/GCA_958296335.1/ensembl/geneset/2024_02/Urtica_urens-GCA_958296335.1-2024_02-cds.fa.gz | 8,781,274 | 956e143f10d1f736bd4d7e3a52edcaf6 | no_md5 |
| 2026-10-23T08:48:27 | UU_ens | Urtica_urens-GCA_958296335.1-2024_02-genes.gff3.gz | https://ftp.ensembl.org/pub/rapid-release/species/Urtica_urens/GCA_958296335.1/ensembl/geneset/2024_02/Urtica_urens-GCA_958296335.1-2024_02-genes.gff3.gz | 4,624,215 | 16420065872feeaa97d02217c4532798 | no_md5 |
| 2026-10-23T08:48:28 | UU_ens | Urtica_urens-GCA_958296335.1-2024_02-pep.fa.gz | https://ftp.ensembl.org/pub/rapid-release/species/Urtica_urens/GCA_958296335.1/ensembl/geneset/2024_02/Urtica_urens-GCA_958296335.1-2024_02-pep.fa.gz | 5,811,135 | fbf2d147752d6264ab87b565fa066e84 | no_md5 |
| 2026-10-23T08:48:28 | UU_ens | md5sum.txt | https://ftp.ensembl.org/pub/rapid-release/species/Urtica_urens/GCA_958296335.1/ensembl/geneset/2024_02/md5sum.txt | 581 | 9fce3db540c95c3581a0ca6bc4d9a07c | no_md5 |
| 2026-10-23T08:48:33 | PJ_ens | Parietaria_judaica-GCA_954871525.1-2024_02-cds.fa.gz | https://ftp.ensembl.org/pub/rapid-release/species/Parietaria_judaica/GCA_954871525.1/ensembl/geneset/2024_02/Parietaria_judaica-GCA_954871525.1-2024_02-cds.fa.gz | 16,995,807 | 07ef10837b9f3adbd66ad74c8b7bc99f | no_md5 |
| 2026-10-23T08:48:35 | PJ_ens | Parietaria_judaica-GCA_954871525.1-2024_02-genes.gff3.gz | https://ftp.ensembl.org/pub/rapid-release/species/Parietaria_judaica/GCA_954871525.1/ensembl/geneset/2024_02/Parietaria_judaica-GCA_954871525.1-2024_02-genes.gff3.gz | 8,749,196 | 43122cd55527a81d188c1f540f4aeab9 | no_md5 |
| 2026-10-23T08:48:37 | PJ_ens | Parietaria_judaica-GCA_954871525.1-2024_02-pep.fa.gz | https://ftp.ensembl.org/pub/rapid-release/species/Parietaria_judaica/GCA_954871525.1/ensembl/geneset/2024_02/Parietaria_judaica-GCA_954871525.1-2024_02-pep.fa.gz | 11,199,888 | a0e61e51b2a5578514d8924fabae2b99 | no_md5 |
| 2026-10-23T08:48:37 | PJ_ens | md5sum.txt | https://ftp.ensembl.org/pub/rapid-release/species/Parietaria_judaica/GCA_954871525.1/ensembl/geneset/2024_02/md5sum.txt | 623 | 62485c40fd0212a50d70256e68e442e6 | no_md5 |

## Supplementary data (fetched 2026-10-23)

| item | source | URL / identifier | local file | use |
|---|---|---|---|---|
| Hirabayashi 2025 supplementary ZIP (Table S8 putative centromere intervals, female H1) | Europe PMC supplementary files for PMC11722821 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11722821/supplementaryFiles (zip containing plants-14-00124-s001.zip → 'Hirabayashi 2024 Supplementary Tables resubmission.xlsx', sheet Table S8) | data/annot/plants2025_supp.zip → data/annot/tableS8_centromere_intervals.csv | H3 centromere points (female) |
| Hirabayashi 2026 supplementary PDF (Supplementary Table 5 centromere intervals, male) | figshare project 230981 (CC BY 4.0) | https://ndownloader.figshare.com/files/61964272 | data/annot/male_supplement.pdf → data/annot/male_tableS5_centromere_intervals.csv (parsed with pypdf; checked by eye) | H3 centromere points (male, held-out) |
| BRAKER3 gene sets (female, male) and EDTA TE annotations | figshare project 230981 (CC BY 4.0) | https://ndownloader.figshare.com/files/51114632, /files/57352571, /files/51114917 (file ids resolved through api.figshare.com project 230981) | data/annot/Figshare_uploads_Gene.tar.gz, Figshare_uploads_Gene_male.tar.gz, Figshare_uploads_TE.tar.gz | inspected; not used (RefSeq annotation preferred) |
| Hirabayashi 2026 Table 2 (maternal/paternal haplotype per chromosome) | Sci Data 13:569 full text (Springer Nature) | https://doi.org/10.1038/s41597-026-06960-7 | transcribed into code/h1_distances.py (MAT = MH1 for chr 1,2,4,7–13; MH2 for chr 3,5,6) | parental labels |
| DToL genome note full text | Wellcome Open Research via CrossRef TDM | https://doi.org/10.12688/wellcomeopenres.23399.1 | articles/10.12688_wellcomeopenres.23399.1.xml | assembly pipeline, k-mer completeness, curation statements |
| BioSample metadata | ENA browser API | https://www.ebi.ac.uk/ena/browser/api/xml/SAMEA7521824 | (read online 2026-10-24) | sex and locality of the DToL plant |

## Software binaries (downloaded 2026-10-22)

| tool | version | URL | size (bytes) | checksum |
|---|---|---|---|---|
| iqtree | iqtree-2.4.0-Windows.zip | https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Windows.zip | 10274929 | 93f488a6cfef96ef664b23fefe5be97c2c80ebcae9cd2031c0cac92b467278be |
| paml | paml-4.10.10-win-x86_64.tar.gz | https://github.com/abacus-gene/paml/releases/download/v4.10.10/paml-4.10.10-win-x86_64.tar.gz | 6168840 | 1de2cd4821391a9b5e32be1699d3afb87c2dc78104275c4dedcb0b3d64a1bdca |
| muscle | muscle-win64.v5.3.exe | https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-win64.v5.3.exe | 1281536 | a86637db99deb3efa1cb6762c4e3160ca33644608270eb167e4214f1f9409711 |
| blast+ | ncbi-blast-2.17.0+-x64-win64.tar.gz | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-win64.tar.gz | 143400333 | md5:dcd973097407a2910061ff4fb51b09fb |

MAFFT not used: mafft.cbrc.jp outside network allowlist; MUSCLE v5.3 substituted.


## Bonn assembly (bonndata doi:10.60507/FK2/1XHSZ3, CC BY 4.0, downloaded 2026-10-23)

| file | dataverse file id | bytes | checksum alg | verified |
|---|---|---|---|---|
| README.txt | 52102 | 4752 | md5 | True |
| Urtica_dioica.v1.cds.fasta.gz | 51392 | 14514285 | md5 | True |
| Urtica_dioica.v1.genome.fasta.gz | 51390 | 274544885 | md5 | True |
| Urtica_dioica.v1.genome.gff.gz | 51394 | 7171721 | md5 | True |
| Urtica_dioica.v1.pep.fasta.gz | 51391 | 9629778 | md5 | True |

Wolff et al. 2026 preprint text (bioRxiv 10.64898/2026.05.15.725508) read online 2026-10-24 for plant origin (Braunschweig; single root; "two haplophases" interpretation).

## Not downloaded
Raw reads (DToL HiFi ERR10662027–ERR10662029; Bonn nanopore PRJEB63450; all 329 Urtica SRA/ENA runs in data/runs.csv) — outside the study's disk/time ceiling; see VALIDATION.md item 5. Ramie (*Boehmeria nivea*) assemblies — no NCBI annotation; not used.
