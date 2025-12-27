suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript deseq2_tximport.R <samples_tsv> <quant_dir> <gtf_gz> <out_tsv>")
}

samples_tsv <- args[[1]]
quant_dir   <- args[[2]]
gtf_gz      <- args[[3]]
out_tsv     <- args[[4]]

# Load sample metadata
samples <- read_tsv(samples_tsv, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id),
         condition = factor(condition))

# Build salmon quant.sf paths
files <- file.path(quant_dir, samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing quant.sf files:\n", paste(missing, collapse = "\n"))
}

# Build tx2gene from GTF
# Keep transcript_id / gene_id, strip version suffixes to match Salmon if needed.
gtf <- read_tsv(gtf_gz, comment = "#", col_names = FALSE, show_col_types = FALSE)
colnames(gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")

tx_rows <- gtf %>% filter(feature == "transcript")

extract_attr <- function(x, key) {
  # extract key "value" from GTF attribute field
  m <- regmatches(x, regexpr(paste0(key, " \"[^\"]+\""), x))
  sub(paste0(key, " \""), "", sub("\"$", "", m))
}

tx2gene <- tx_rows %>%
  transmute(
    tx = extract_attr(attribute, "transcript_id"),
    gene = extract_attr(attribute, "gene_id")
  ) %>%
  filter(!is.na(tx), !is.na(gene)) %>%
  mutate(
    tx = sub("\\..*$", "", tx),
    gene = sub("\\..*$", "", gene)
  ) %>%
  distinct()

# Import Salmon
# countsFromAbundance helps DESeq2 by using length-scaled counts. :contentReference[oaicite:2]{index=2}
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE,
  ignoreAfterBar = TRUE
)

# DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = as.data.frame(samples), design = ~ condition)
dds <- DESeq(dds)

# Treated vs control (make sure factor levels are correct)
# If your levels are control/treated, this yields treated_vs_control.
res <- results(dds, contrast = c("condition", "treated", "control"))

out <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(padj)

write_tsv(out, out_tsv)
