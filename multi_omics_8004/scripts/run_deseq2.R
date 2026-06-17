#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis
# Usage: Rscript run_deseq2.R <input_csv> <output_dir> [options_json]
# Input: CSV with genes x samples (first col = gene names)
# Output: DE results CSV, volcano plot PNG, normalized counts, MA plot

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat('Usage: Rscript run_deseq2.R <input_csv> <output_dir> [options_json]\n')
  quit(status = 1)
}

input_file <- args[1]
output_dir <- args[2]
options_json <- if (length(args) >= 3) jsonlite::fromJSON(args[3]) else list()

# Parameters with defaults
pval_cutoff <- if (!is.null(options_json$pvalue_cutoff)) options_json$pvalue_cutoff else 0.05
log2fc_cutoff <- if (!is.null(options_json$log2fc_cutoff)) options_json$log2fc_cutoff else 1.0
control_prefix <- if (!is.null(options_json$control_prefix)) options_json$control_prefix else 'ctrl'
treat_prefix <- if (!is.null(options_json$treat_prefix)) options_json$treat_prefix else 'treat'

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf('[INFO] Reading input: %s\n', input_file))
df <- read.csv(input_file, check.names = FALSE, row.names = 1)

# Find gene column if rownames not set
if ('gene' %in% tolower(colnames(df))) {
  gene_col <- which(tolower(colnames(df)) == 'gene')
  rownames(df) <- df[, gene_col]
  df <- df[, -gene_col]
}

# Keep only numeric columns
num_cols <- sapply(df, is.numeric)
df <- df[, num_cols, drop = FALSE]
df <- as.data.frame(lapply(df, as.integer))
rownames(df) <- make.unique(rownames(df))

cat(sprintf('[INFO] Data: %d genes x %d samples\n', nrow(df), ncol(df)))
cat(sprintf('[INFO] Samples: %s\n', paste(colnames(df), collapse = ', ')))

# Determine groups from column names
# Strategy: first half = control, second half = treatment (or use prefix matching)
conditions <- rep('treatment', ncol(df))
ctrl_match <- grep(control_prefix, colnames(df), ignore.case = TRUE)
treat_match <- grep(treat_prefix, colnames(df), ignore.case = TRUE)

if (length(ctrl_match) > 0 && length(treat_match) > 0) {
  conditions[ctrl_match] <- 'control'
  conditions[treat_match] <- 'treatment'
} else {
  # Default: first half control, second half treatment
  mid <- ncol(df) %/% 2
  conditions[1:mid] <- 'control'
  if (mid < ncol(df)) conditions[(mid+1):ncol(df)] <- 'treatment'
}

coldata <- data.frame(
  condition = factor(conditions, levels = c('control', 'treatment')),
  row.names = colnames(df)
)

cat(sprintf('[INFO] Groups: %d control, %d treatment\n',
    sum(conditions == 'control'), sum(conditions == 'treatment')))

# DESeq2 analysis
cat('[INFO] Running DESeq2...\n')
dds <- DESeqDataSetFromMatrix(countData = df, colData = coldata, design = ~ condition)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= min(2, ncol(dds) %/% 2)
if (sum(keep) < nrow(dds)) {
  cat(sprintf('[INFO] Filtered %d low-count genes (%d remaining)\n', nrow(dds) - sum(keep), sum(keep)))
}
dds <- dds[keep, ]

dds <- DESeq(dds)

# Results
res <- results(dds, alpha = pval_cutoff, contrast = c('condition', 'treatment', 'control'))
# lfcShrink skipped - using raw results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# Add significance column
res_df$significant <- ifelse(is.na(res_df$padj), 0,
                      ifelse(abs(res_df$log2FoldChange) >= log2fc_cutoff & res_df$padj < pval_cutoff, 1, 0))

n_total <- nrow(res_df)
n_sig <- sum(res_df$significant == 1)
n_up <- sum(res_df$significant == 1 & res_df$log2FoldChange > 0)
n_down <- sum(res_df$significant == 1 & res_df$log2FoldChange < 0)

cat(sprintf('[RESULT] %d genes analyzed, %d significant (%d up, %d down)\n', n_total, n_sig, n_up, n_down))

# Save DE results
write.csv(res_df, file.path(output_dir, 'deseq2_results.csv'), row.names = FALSE)
cat(sprintf('[INFO] Results saved to %s/deseq2_results.csv\n', output_dir))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_df <- as.data.frame(norm_counts)
norm_df$gene <- rownames(norm_df)
write.csv(norm_df, file.path(output_dir, 'normalized_counts.csv'), row.names = FALSE)

# --- Volcano Plot (EnhancedVolcano) ---
cat('[INFO] Generating volcano plot...\n')
png(file.path(output_dir, 'volcano_plot.png'), width = 3000, height = 2400, res = 300)
p <- EnhancedVolcano(res_df,
    lab = res_df$gene,
    x = 'log2FoldChange',
    y = 'padj',
    title = 'DESeq2 Differential Expression',
    subtitle = sprintf('control vs treatment (padj < %s, |log2FC| > %s)', pval_cutoff, log2fc_cutoff),
    pCutoff = pval_cutoff,
    FCcutoff = log2fc_cutoff,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.8,
    legendPosition = 'right',
    border = 'full',
    borderWidth = 1.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE)
print(p)
dev.off()

# --- MA Plot ---
cat('[INFO] Generating MA plot...\n')
png(file.path(output_dir, 'ma_plot.png'), width = 3000, height = 2400, res = 300)
plotMA(res, main = 'DESeq2 MA Plot', ylim = c(-5, 5))
dev.off()

# --- Top gene heatmap ---
cat('[INFO] Generating heatmap...\n')
top_genes <- head(rownames(res[order(res$padj), ]), min(30, nrow(res)))
if (length(top_genes) >= 2) {
  mat <- norm_counts[top_genes, ]
  mat <- mat[complete.cases(mat), ]
  if (nrow(mat) >= 2) {
    mat_scaled <- t(scale(t(mat)))
    png(file.path(output_dir, 'heatmap.png'), width = 2400, height = 3000, res = 300)
    pheatmap(mat_scaled,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             main = sprintf('Top %d DE Genes', nrow(mat_scaled)),
             annotation_col = coldata,
             fontsize_row = 8)
    dev.off()
  }
}

# Summary JSON
summary <- list(
  total_genes = n_total,
  significant = n_sig,
  up_regulated = n_up,
  down_regulated = n_down,
  pval_cutoff = pval_cutoff,
  log2fc_cutoff = log2fc_cutoff,
  n_control = sum(conditions == 'control'),
  n_treatment = sum(conditions == 'treatment'),
  method = 'DESeq2',
  output_files = c('deseq2_results.csv', 'normalized_counts.csv',
                   'volcano_plot.png', 'ma_plot.png', 'heatmap.png')
)
write(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
      file.path(output_dir, 'summary.json'))

cat('[DONE] DESeq2 analysis complete!\n')
