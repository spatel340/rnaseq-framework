# RNA-seq End-to-End Framework (FASTQ → QC → DE → Enrichment → Report)

This repository is a reproducible RNA-seq analysis workflow built with **Snakemake**.
Milestone 1 implements: **SRA download → FASTQ generation → FastQC → MultiQC**.

## Quickstart
1) Edit `data/metadata/samples.tsv` to include real SRA run accessions (SRR…).
2) Run QC pipeline:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 4

q
