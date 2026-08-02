# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/utils.R
# ---------------------------------------------------------------------------
# Shared helper functions used across all pipeline stages.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
})

# Create a directory if it doesn't already exist.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

# Save a ggplot with consistent size/resolution. bg = "white" avoids
# transparent panels (e.g. cnetplot, theme_void) rendering as solid black
# in some viewers.
save_plot <- function(plot, filename, width = 7, height = 6, dpi = 300) {
  ensure_dir(dirname(filename))
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
}

# Read a cached .rds if present, otherwise evaluate `expr` and cache it.
# Lets any step be re-run without repeating expensive work.
cache_step <- function(path, expr) {
  if (file.exists(path)) return(readRDS(path))
  result <- expr
  ensure_dir(dirname(path))
  saveRDS(result, path)
  result
}

# Timestamped progress message, printed to console.
log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# Boxplot of normalized expression for one gene, split by sample condition.
plot_gene_of_interest <- function(vsd, gene_info, gene_symbol, col_data) {
  gene_id <- gene_info$gene_id[gene_info$gene_name == gene_symbol][1]
  if (is.na(gene_id)) stop("Gene symbol not found: ", gene_symbol)

  expr <- SummarizedExperiment::assay(vsd)[gene_id, ]
  df <- data.frame(expression = expr, condition = col_data$condition)

  ggplot(df, aes(condition, expression, fill = condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4) +
    labs(title = gene_symbol, y = "VST-normalized expression", x = NULL) +
    theme_bw() +
    theme(legend.position = "none")
}

# Pathway/term names from MSigDB, GO, KEGG, and Reactome share a long
# collection prefix (e.g. "GOBP_") and use underscores instead of spaces.
# Strips the shared prefix, title-cases, and wraps to a fixed width so
# labels stay legible in plots. Shared by GSEA, GO, and Reactome plotting.
clean_pathway_labels <- function(x, wrap_width = 40) {
  strip_common_prefix <- function(labels) {
    if (length(labels) <= 1) return(labels)
    shortest <- min(nchar(labels))
    prefix_len <- 0
    for (i in seq_len(shortest)) {
      if (length(unique(substr(labels, i, i))) > 1) break
      prefix_len <- i
    }
    if (prefix_len >= 5) substring(labels, prefix_len + 1) else labels
  }

  x <- strip_common_prefix(x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  str_wrap(x, width = wrap_width)
}

# Opens a file in the system's default browser/viewer. Used to auto-launch
# the HTML report after rendering. Never fails the pipeline -- if no
# display is available (e.g. a headless server), it logs and continues.
launch_file <- function(path) {
  tryCatch({
    utils::browseURL(normalizePath(path))
  }, error = function(e) {
    log_step(paste("Could not auto-launch report --", conditionMessage(e)))
  })
}
