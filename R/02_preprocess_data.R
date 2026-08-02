# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/02_preprocess_data.R
# ---------------------------------------------------------------------------
# Turns the raw SummarizedExperiment into the plain objects every
# downstream script needs: a filtered counts matrix, sample metadata, a
# gene-ID/symbol lookup, and processed clinical data (with optional
# covariate collapsing via config.R).
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# --- Pre-filtering ---------------------------------------------------------

se <- readRDS(file.path(DATA_DIR, "raw_se.rds"))
# head(se)    # uncomment to explore

# Gene ID <-> symbol lookup, cached before `se` goes out of scope so no
# downstream script depends on the (large) SummarizedExperiment staying
# in memory.
gene_info <- cache_step(file.path(DATA_DIR, "gene_info.rds"), {
  as.data.frame(rowData(se)) %>%
    select(gene_id, gene_name) %>%
    distinct()
})

counts <- assay(se, "unstranded")

# Sample condition (Tumor / Normal) inferred from the TCGA barcode's sample
# type code -- codes 10/11 (Blood/Solid Tissue Normal) are Normal, all
# else Tumor. Same logic for any TCGA project.
sample_type_code <- substr(colnames(counts), 14, 15)
condition <- ifelse(sample_type_code %in% c("10", "11"), "Normal", "Tumor")

col_data <- data.frame(
  barcode   = colnames(counts),
  condition = factor(condition),
  row.names = colnames(counts)
)

# Drop lowly-expressed genes. Threshold is relative to sample count rather
# than a fixed number, so it scales to datasets of any size.
keep <- rowSums(counts >= 10) >= 0.2 * ncol(counts)
counts <- counts[keep, ]

n_tumor  <- sum(col_data$condition == "Tumor")
n_normal <- sum(col_data$condition == "Normal")

log_step(sprintf(
  "Preprocessed: %d genes retained, %d Tumor / %d Normal samples.",
  nrow(counts), n_tumor, n_normal
))

# Some TCGA projects have very few or zero matched Normal samples, which
# makes the Tumor-vs-Normal DESeq2 contrast statistically fragile. This
# won't stop the pipeline, but it's worth knowing before trusting the DEG
# list downstream.
if (n_normal < MIN_NORMAL_SAMPLES) {
  log_step(sprintf(
    "WARNING: only %d Normal samples (< MIN_NORMAL_SAMPLES = %d). Tumor-vs-Normal DESeq2 results may be underpowered/unstable for this project.",
    n_normal, MIN_NORMAL_SAMPLES
  ))
}

saveRDS(counts, file.path(DATA_DIR, "counts_filtered.rds"))
saveRDS(col_data, file.path(DATA_DIR, "col_data.rds"))

# --- Clinical data processing ------------------------------------------------
# Collapses sparse factor levels (e.g. AJCC sub-stages) per config.R's
# COVARIATE_COLLAPSE_MAP, so downstream Cox models don't hit separation.
# A covariate with no map entry is left as-is, so this is a no-op for
# projects that don't define one.
clinical <- readRDS(file.path(DATA_DIR, "clinical.rds"))
collapse_map <- if (exists("COVARIATE_COLLAPSE_MAP")) COVARIATE_COLLAPSE_MAP else list()

for (cov in intersect(names(collapse_map), colnames(clinical))) {
  # recode() leaves values not listed in the map unchanged.
  clinical[[cov]] <- dplyr::recode(as.character(clinical[[cov]]), !!!collapse_map[[cov]])
}

saveRDS(clinical, file.path(DATA_DIR, "clinical_processed.rds"))

log_step(sprintf(
  "Clinical data processed: %d patients, %d covariate(s) collapsed.",
  nrow(clinical), length(intersect(names(collapse_map), colnames(clinical)))
))
log_step("02_preprocess_data.R finished.")
