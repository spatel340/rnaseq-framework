# RNA-seq End-to-End Analysis Framework (Snakemake)

## Problem
Quantify transcriptional changes between **treated vs control** samples using bulk RNA-seq and generate interpretable biological conclusions.

## Data
Public RNA-seq dataset from NCBI SRA (run accessions listed in `data/metadata/samples.tsv`).
This workflow downloads raw data programmatically (no manual downloads required).

## Methods
Pipeline stages (all automated with Snakemake + conda envs):
1. **Data acquisition**: `prefetch` (SRA) → `fasterq-dump` (FASTQ)
2. **Quality control**: FastQC + MultiQC
3. **Quantification**: Salmon (transcript-level) + gene summarization via tximport
4. **Differential expression**: DESeq2 (`~ condition`)
5. **Interpretation**: GO Biological Process enrichment
6. **Reporting**: Volcano + GO dotplot figures

## Results (high level)
Artifacts are produced under `results/`:
- `results/qc/multiqc_report.html`
- `results/de/deseq2_results.tsv`
- `results/figures/volcano.pdf`
- `results/figures/go_bp_dotplot.pdf`
- `results/enrichment/go_bp.tsv`

**Key takeaways**
- [Add 2–3 bullets after you inspect volcano + GO plot]

## How to reproduce
### Prerequisites
- Linux/WSL
- Conda/Miniforge

### Run
```bash
snakemake -s workflow/Snakefile --use-conda --cores 2 --jobs 1
