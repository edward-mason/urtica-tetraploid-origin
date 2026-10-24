# Tetraploid origin of *Urtica dioica* — preregistered assembly-based test

Companion repository for: Mason E. (2026) *No second subgenome more divergent than diploid alleles in the deposited genome of a tetraploid stinging nettle (Urtica dioica): a preregistered genome-wide test against haplotype-resolved diploid assemblies.* bioRxiv preprint.

- `paper/` — manuscript (Markdown and submission PDF), figures, claim ledger, AI disclosure, independent reviewer reports and responses
- `notes/` — preregistration (v1/v2 with SHA-256), hypothesis list, laboratory notebook, literature map
- `code/` — Snakemake workflow and Python scripts (see `code/README.md`); BLAST+ 2.17.0 required
- `data/` — provenance (every download with MD5), inventories of public assemblies/annotations/reads, centromere tables derived from published supplements. Genome FASTA files are not included (download with `code/download.py`).
- `results/` — per-gene distance tables, H1/H2/H3 result files, synteny blocks, sensitivity tables. Raw BLAST tables (~5 GB) are not included; they regenerate from the workflow.

Licence: code MIT; text, figures and tables CC BY 4.0. All input data are public (NCBI, Ensembl, bonndata, figshare, Europe PMC); see `data/PROVENANCE.md`.
