# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/04_functional_enrichment.R
# ---------------------------------------------------------------------------
# Over-representation analysis (GO + KEGG + Reactome) on the significant
# DEGs, run through one function so all three share the same Entrez ID
# mapping, background/universe, and significance threshold.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(dplyr)
})

res <- readRDS(file.path(DATA_DIR, "deseq2_res.rds"))

run_ora <- function(gene_symbols, universe_symbols, direction_label) {
  id_map       <- clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = ORGANISM_DB)
  universe_map <- clusterProfiler::bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = ORGANISM_DB)

  log_step(sprintf(
    "%s: %d/%d genes, %d/%d universe genes mapped to Entrez.",
    direction_label, nrow(id_map), length(gene_symbols), nrow(universe_map), length(universe_symbols)
  ))

  # pvalueCutoff is tied to config.R's PVAL_CUTOFF so the enrichment
  # significance threshold moves with the rest of the pipeline instead of
  # silently staying at clusterProfiler's own default (0.05) if PVAL_CUTOFF
  # is changed. qvalueCutoff = 1 leaves FDR q-value filtering off; padj is
  # still reported in the output table for the reader to filter on.
  go <- tryCatch({
    clusterProfiler::enrichGO(
      gene = id_map$ENTREZID, universe = universe_map$ENTREZID,
      OrgDb = ORGANISM_DB, ont = "BP", pAdjustMethod = "BH", readable = TRUE,
      pvalueCutoff = PVAL_CUTOFF, qvalueCutoff = 1
    )
  }, error = function(e) {
    log_step(paste("GO enrichment failed (", direction_label, ") --", conditionMessage(e)))
    NULL
  })

  # enrichKEGG() fetches live from rest.kegg.jp -- the default 60s download
  # timeout can be too short for a slow response. Raised here and restored
  # after, with one retry at a longer timeout if the first attempt fails.
  fetch_kegg <- function(timeout_secs) {
    old_timeout <- getOption("timeout")
    options(timeout = timeout_secs)
    on.exit(options(timeout = old_timeout), add = TRUE)
    kegg_raw <- clusterProfiler::enrichKEGG(
      gene = id_map$ENTREZID, universe = universe_map$ENTREZID,
      organism = KEGG_ORGANISM, pAdjustMethod = "BH",
      pvalueCutoff = PVAL_CUTOFF, qvalueCutoff = 1
    )
    clusterProfiler::setReadable(kegg_raw, OrgDb = ORGANISM_DB, keyType = "ENTREZID")
  }

  kegg <- tryCatch(fetch_kegg(120), error = function(e) {
    log_step(paste("KEGG fetch failed at 120s, retrying at 300s --", conditionMessage(e)))
    tryCatch(fetch_kegg(300), error = function(e2) {
      log_step(paste("KEGG enrichment failed (", direction_label, ") --", conditionMessage(e2)))
      NULL
    })
  })

  # Namespace-qualified: enrichPathway can otherwise fail to resolve if
  # another loaded package's search-path position shadows it.
  reactome <- tryCatch({
    ReactomePA::enrichPathway(
      gene = id_map$ENTREZID, universe = universe_map$ENTREZID,
      organism = REACTOME_ORGANISM, pAdjustMethod = "BH", readable = TRUE,
      pvalueCutoff = PVAL_CUTOFF, qvalueCutoff = 1
    )
  }, error = function(e) {
    log_step(paste("Reactome enrichment failed (", direction_label, ") --", conditionMessage(e)))
    NULL
  })

  list(go = go, kegg = kegg, reactome = reactome, direction = direction_label)
}

save_ora_outputs <- function(ora) {
  tag <- tolower(ora$direction)

  if (!is.null(ora$go)) {
    write.csv(as.data.frame(ora$go), file.path(RESULTS_DIR, paste0("go_bp_", tag, ".csv")), row.names = FALSE)
  }
  if (!is.null(ora$kegg)) {
    write.csv(as.data.frame(ora$kegg), file.path(RESULTS_DIR, paste0("kegg_", tag, ".csv")), row.names = FALSE)
  }
  if (!is.null(ora$reactome)) {
    write.csv(as.data.frame(ora$reactome), file.path(RESULTS_DIR, paste0("reactome_", tag, ".csv")), row.names = FALSE)
  }

  # dotplot()'s label_format wraps long term names; height scales with the
  # number of terms shown so wrapped labels don't collide vertically.
  if (!is.null(ora$go)) {
    n_go <- min(15, nrow(as.data.frame(ora$go)))
    if (n_go > 0) {
      go_plot <- dotplot(ora$go, showCategory = n_go, label_format = 45) +
        ggplot2::labs(title = paste("GO Biological Process -", ora$direction), subtitle = PROJECT) +
        ggplot2::theme(plot.title.position = "plot", plot.margin = ggplot2::margin(12, 20, 10, 10))
      save_plot(go_plot, file.path(RESULTS_DIR, paste0("go_bp_", tag, ".png")),
                width = 10, height = max(6, 0.45 * n_go + 2.5))
    }
  }

  if (!is.null(ora$kegg)) {
    save_plot(
      cnetplot(ora$kegg, showCategory = 8) +
        ggplot2::labs(title = paste("KEGG Pathways -", ora$direction), subtitle = PROJECT) +
        ggplot2::theme(plot.title.position = "plot", plot.margin = ggplot2::margin(12, 20, 10, 10)),
      file.path(RESULTS_DIR, paste0("kegg_", tag, ".png")), width = 11, height = 9
    )
  }

  if (!is.null(ora$reactome)) {
    n_reactome <- min(15, nrow(as.data.frame(ora$reactome)))
    if (n_reactome > 0) {
      reactome_plot <- dotplot(ora$reactome, showCategory = n_reactome, label_format = 45) +
        ggplot2::labs(title = paste("Reactome Pathways -", ora$direction), subtitle = PROJECT) +
        ggplot2::theme(plot.title.position = "plot", plot.margin = ggplot2::margin(12, 20, 10, 10))
      save_plot(reactome_plot, file.path(RESULTS_DIR, paste0("reactome_", tag, ".png")),
                width = 10, height = max(6, 0.45 * n_reactome + 2.5))
    }
  }

  log_step(sprintf(
    "%s enrichment saved: GO %s, KEGG %s, Reactome %s.",
    ora$direction,
    if (is.null(ora$go)) "skipped" else "ok",
    if (is.null(ora$kegg)) "skipped" else "ok",
    if (is.null(ora$reactome)) "skipped" else "ok"
  ))
}

up_genes   <- res %>% filter(padj < PVAL_CUTOFF, log2FoldChange >=  LFC_CUTOFF) %>% pull(gene_name)
down_genes <- res %>% filter(padj < PVAL_CUTOFF, log2FoldChange <= -LFC_CUTOFF) %>% pull(gene_name)
background <- res %>% pull(gene_name)

save_ora_outputs(run_ora(up_genes, background, "Upregulated"))
save_ora_outputs(run_ora(down_genes, background, "Downregulated"))

log_step("04_functional_enrichment.R finished.")
