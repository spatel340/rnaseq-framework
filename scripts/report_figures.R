suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript report_figures.R <deseq2_results_tsv>")
res_path <- args[[1]]

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

res <- read_tsv(res_path, show_col_types = FALSE)

# Basic filtering
res2 <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(sig = padj < 0.05)

# Volcano
pdf("results/figures/volcano.pdf", width = 7, height = 6)
print(
  EnhancedVolcano(
    res2,
    lab = res2$gene_id,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 1.5,
    labSize = 2.0,
    title = "DESeq2: treated vs control",
    subtitle = "Volcano plot (padj < 0.05, |log2FC| > 1)"
  )
)
dev.off()

# Enrichment: map gene_id -> Entrez
# If your gene_id is Ensembl gene IDs, map to Entrez.
entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = res2$gene_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res2$ENTREZID <- unname(entrez)

sig_entrez <- res2 %>%
  filter(sig, !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

if (length(sig_entrez) < 10) {
  message("Too few significant genes mapped to Entrez for enrichment; skipping GO.")
  writeLines("Too few significant genes for GO enrichment.", "results/enrichment/README.txt")
  quit(status = 0)
}

ego <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

write_tsv(as_tibble(ego@result), "results/enrichment/go_bp.tsv")

pdf("results/figures/go_bp_dotplot.pdf", width = 8, height = 6)
print(dotplot(ego, showCategory = 15) + ggtitle("GO BP enrichment (sig DE genes)"))
dev.off()
