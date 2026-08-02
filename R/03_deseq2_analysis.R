# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/03_deseq2_analysis.R
# ---------------------------------------------------------------------------
# Differential expression with DESeq2, plus the standard accompanying
# QC/results plots: PCA, dispersion, volcano, MA plot, and a top-DEG
# heatmap.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(ashr)
})

ensure_dir(RESULTS_DIR)

counts    <- readRDS(file.path(DATA_DIR, "counts_filtered.rds"))
col_data  <- readRDS(file.path(DATA_DIR, "col_data.rds"))
gene_info <- readRDS(file.path(DATA_DIR, "gene_info.rds"))

# --- Run DESeq ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
dds <- DESeq(dds)

# --- Results -----------------------------------------------------------------
res_raw <- results(dds, contrast = c("condition", "Tumor", "Normal"))

res <- res_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by = "gene_id") %>%
  filter(!is.na(padj))

sig_genes <- res %>%
  filter(padj < PVAL_CUTOFF, abs(log2FoldChange) >= LFC_CUTOFF)

log_step(sprintf(
  "DESeq2 complete: %d significant genes (padj < %s, |LFC| >= %s).",
  nrow(sig_genes), PVAL_CUTOFF, LFC_CUTOFF
))

write.csv(res, file.path(RESULTS_DIR, "deseq2_results.csv"), row.names = FALSE)

# Variance-stabilised counts -- reused by the QC plots below, GSEA ranking,
# and the survival analysis script.
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(DATA_DIR, "vsd.rds"))
saveRDS(res, file.path(DATA_DIR, "deseq2_res.rds"))

# --- Visualization -----------------------------------------------------------

# --- PCA ---
tryCatch({
  pca_plot <- plotPCA(vsd, intgroup = "condition") + theme_bw()
  save_plot(pca_plot, file.path(RESULTS_DIR, "qc_pca.png"))
  log_step("PCA plot generated.")
}, error = function(e) {
  log_step(paste("PCA plot failed --", conditionMessage(e)))
})

# --- Dispersion plot (always; independent of Tumor/Normal group sizes) ---
tryCatch({
  png(file.path(RESULTS_DIR, "qc_dispersion.png"), width = 7, height = 6, units = "in", res = 300)
  plotDispEsts(dds, main = paste(PROJECT, "- Dispersion Estimates"))
  dev.off()
  log_step("Dispersion plot generated.")
}, error = function(e) {
  if (dev.cur() != 1) dev.off()
  log_step(paste("Dispersion plot failed --", conditionMessage(e)))
})

# --- Poisson distance heatmap (QC fallback for low/no Normal samples) ---
# Sample-to-sample clustering on raw counts, independent of group label --
# catches outliers/batch structure when there's too little Normal group
# for the PCA/condition split above to be informative.
n_normal <- sum(col_data$condition == "Normal")
if (n_normal < MIN_NORMAL_SAMPLES) {
  tryCatch({
    if (!requireNamespace("PoiClaClu", quietly = TRUE)) {
      stop("PoiClaClu not installed (install.packages('PoiClaClu'))")
    }
    poisd     <- PoiClaClu::PoissonDistance(t(counts))
    poisd_mat <- as.matrix(poisd$dd)
    dimnames(poisd_mat) <- list(colnames(counts), colnames(counts))

    pheatmap(
      poisd_mat,
      clustering_distance_rows = poisd$dd,
      clustering_distance_cols = poisd$dd,
      annotation_col = col_data["condition"],
      show_rownames  = FALSE, show_colnames = FALSE,
      main           = paste(PROJECT, "- Poisson Sample Distance (QC fallback)"),
      filename       = file.path(RESULTS_DIR, "qc_poisson_distance.png")
    )
    log_step(sprintf("Poisson distance heatmap generated (only %d Normal samples).", n_normal))
  }, error = function(e) {
    log_step(paste("Poisson distance heatmap skipped --", conditionMessage(e)))
  })
}

# --- Volcano ---
tryCatch({
  top_labels <- res %>% arrange(padj) %>% slice_head(n = 10) %>% pull(gene_name)
  if (nzchar(GENE_OF_INTEREST)) top_labels <- union(top_labels, GENE_OF_INTEREST)

  volcano <- EnhancedVolcano(
    res,
    lab       = res$gene_name,
    x         = "log2FoldChange",
    y         = "padj",
    selectLab = top_labels,
    pCutoff   = PVAL_CUTOFF,
    FCcutoff  = LFC_CUTOFF,
    title     = paste(PROJECT, "- Tumor vs Normal")
  )
  save_plot(volcano, file.path(RESULTS_DIR, "deg_volcano.png"), width = 8, height = 9)
  log_step("Volcano plot generated.")
}, error = function(e) {
  log_step(paste("Volcano plot failed --", conditionMessage(e)))
})

# --- MA plot ---
# log2FoldChange vs. mean normalized expression -- checks for
# expression-level-dependent bias and shows how many genes were called
# significant at each expression level. Uses shrunken LFC (config.R:
# LFC_SHRINK_METHOD) if set -- shrinkage only affects this plot;
# `res`/deseq2_results.csv stay unshrunken since GSEA ranks by `stat`,
# which shrinkage doesn't provide.
#
# plotMA() requires a DESeqResults/DESeqDataSet object, not a data.frame --
# both the shrunken result and the unshrunken fallback below use res_raw
# (the S4 object), never the flattened data.frame `res`.
ma_input <- res_raw
ma_shrunk <- FALSE

if (exists("LFC_SHRINK_METHOD") && !identical(LFC_SHRINK_METHOD, "none") && nzchar(LFC_SHRINK_METHOD)) {
  ma_input <- tryCatch({
    shrunk <- lfcShrink(dds, contrast = c("condition", "Tumor", "Normal"), res = res_raw, type = LFC_SHRINK_METHOD)
    ma_shrunk <- TRUE
    shrunk
  }, error = function(e) {
    log_step(paste("LFC shrinkage failed, using unshrunken results for MA plot --", conditionMessage(e)))
    res_raw
  })
}

tryCatch({
  png(file.path(RESULTS_DIR, "qc_ma_plot.png"), width = 7, height = 6, units = "in", res = 300)
  ma_title <- paste(PROJECT, "- MA Plot (Tumor vs Normal)", if (ma_shrunk) sprintf("[%s-shrunken LFC]", LFC_SHRINK_METHOD) else "")
  DESeq2::plotMA(ma_input, alpha = PVAL_CUTOFF, main = ma_title)
  dev.off()
  log_step(paste("MA plot generated.", if (ma_shrunk) sprintf("(LFC shrunk via %s)", LFC_SHRINK_METHOD) else ""))
}, error = function(e) {
  if (dev.cur() != 1) dev.off()
  log_step(paste("MA plot failed --", conditionMessage(e)))
})

# --- Top-DEG heatmap ---
tryCatch({
  top_deg_ids <- sig_genes %>% arrange(padj) %>% slice_head(n = 50) %>% pull(gene_id)
  mat <- assay(vsd)[top_deg_ids, ]
  mat <- t(scale(t(mat)))  # z-score by gene
  rownames(mat) <- gene_info$gene_name[match(rownames(mat), gene_info$gene_id)]

  pheatmap(
    mat,
    annotation_col = col_data["condition"],
    show_colnames  = FALSE,
    filename       = file.path(RESULTS_DIR, "deg_top50_heatmap.png")
  )
  log_step(sprintf("Top-DEG heatmap generated (%d genes).", length(top_deg_ids)))
}, error = function(e) {
  log_step(paste("Top-DEG heatmap failed --", conditionMessage(e)))
})

# --- Gene-of-interest expression plot ---
tryCatch({
  goi_plot <- plot_gene_of_interest(vsd, gene_info, GENE_OF_INTEREST, col_data)
  save_plot(goi_plot, file.path(RESULTS_DIR, paste0("expression_", GENE_OF_INTEREST, ".png")), width = 5, height = 6)
  log_step(paste("Gene-of-interest expression plot generated for", GENE_OF_INTEREST))
}, error = function(e) {
  log_step(paste("Gene-of-interest plot failed --", conditionMessage(e)))
})

log_step("03_deseq2_analysis.R finished.")
