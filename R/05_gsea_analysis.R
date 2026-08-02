# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/05_gsea_analysis.R
# ---------------------------------------------------------------------------
# Gene Set Enrichment Analysis on the full ranked gene list, using MSigDB
# gene sets via msigdbr so the reference collection can be swapped in
# config.R without touching this script.
#
# Uses msigdbr's collection/subcollection API. Run
# msigdbr::msigdbr_collections() to see all valid combinations for your
# installed msigdbr version.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(ggplot2)
})

res <- readRDS(file.path(DATA_DIR, "deseq2_res.rds"))

# Rank genes by DESeq2's test statistic -- more robust to noisy fold
# changes from low-count genes than ranking by log2FoldChange alone.
ranked_genes <- res %>%
  filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  arrange(desc(stat))

gene_ranks <- setNames(ranked_genes$stat, ranked_genes$gene_name)

gene_sets <- msigdbr(
  species       = "Homo sapiens",
  collection    = GSEA_COLLECTION,
  subcollection = if (nzchar(GSEA_SUBCOLLECTION)) GSEA_SUBCOLLECTION else NULL
) %>%
  split(x = .$gene_symbol, f = .$gs_name)

gsea_res <- fgsea(pathways = gene_sets, stats = gene_ranks, minSize = 15, maxSize = 500) %>%
  arrange(padj) %>%
  # Adds a readable label alongside the raw `pathway` column (e.g. "Origin
  # Unwinding And Elongation" instead of "KEGG_MEDICUS_REFERENCE_ORIGIN_...")
  # for display; `pathway` itself stays untouched for downstream analysis.
  mutate(pathway_label = clean_pathway_labels(pathway, wrap_width = 200))

write.csv(gsea_res %>% dplyr::select(-leadingEdge), file.path(RESULTS_DIR, "gsea_results.csv"), row.names = FALSE)

top_pathways <- gsea_res %>% filter(padj < PVAL_CUTOFF) %>% slice_head(n = 15)

if (nrow(top_pathways) > 0) {
  geneset_label <- if (nzchar(GSEA_SUBCOLLECTION)) GSEA_SUBCOLLECTION else GSEA_COLLECTION

  # Wrap y-axis pathway names (utils.R) -- collection names are long, and
  # without wrapping they get truncated or squeeze the plot panel down.
  top_pathways <- top_pathways %>%
    mutate(pathway_label = clean_pathway_labels(pathway, wrap_width = 42))

  # Wrap title/subtitle the same way; legend on the bottom (horizontal)
  # rather than the right so it doesn't eat into the title's width.
  plot_title    <- paste(PROJECT, "- Gene Set Enrichment Analysis")
  plot_subtitle <- paste("Gene set:", geneset_label)

  gsea_plot <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_label, NES), fill = padj)) +
    geom_col() +
    scale_fill_gradient(
      low = "firebrick", high = "steelblue", name = "Adjusted p-value",
      labels = scales::label_scientific()
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Normalized Enrichment Score", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.margin  = margin(t = 12, r = 20, b = 10, l = 10),
      axis.text.y  = element_text(size = 9, lineheight = 0.9),
      legend.position = "bottom",
      legend.key.width = unit(1.4, "cm")
    )

  # Height scales with the number of bars; width is fixed since labels are
  # already wrapped to a known character width.
  plot_height <- max(6, 0.6 * nrow(top_pathways) + 2.5)
  save_plot(gsea_plot, file.path(RESULTS_DIR, "gsea_top_pathways.png"), width = 11, height = plot_height)
}

log_step(sprintf(
  "GSEA complete: %d pathways tested, %d significant at padj < %s.",
  nrow(gsea_res), nrow(top_pathways), PVAL_CUTOFF
))
log_step("05_gsea_analysis.R finished.")
