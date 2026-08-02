# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/01_download_data.R
# ---------------------------------------------------------------------------
# Download (or load from cache) TCGA RNA-seq and clinical data via
# TCGAbiolinks. Change PROJECT in config.R to point at a different
# project -- nothing here needs to change.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

ensure_dir(DATA_DIR)

raw_download_dir <- file.path(DATA_DIR, PROJECT)

se_path       <- file.path(DATA_DIR, "raw_se.rds")
clinical_path <- file.path(DATA_DIR, "clinical.rds")

# --- RNA-seq counts ---------------------------------------------------------
if (!file.exists(se_path)) {
  log_step(paste("Querying GDC for", PROJECT, "RNA-seq data..."))

  rna_query <- GDCquery(
    project       = PROJECT,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = SAMPLE_TYPES,
    access        = "open"
  )

  GDCdownload(rna_query, method = "api", files.per.chunk = 6, directory = DATA_DIR)
  se <- GDCprepare(rna_query, directory = DATA_DIR, summarizedExperiment = TRUE)

  saveRDS(se, se_path)

  # Raw GDC files are only needed to build `se`; once it's cached as
  # raw_se.rds every downstream script reads that file, never the raw
  # download, so it's safe to remove here to save disk space.
  if (dir.exists(raw_download_dir)) {
    unlink(raw_download_dir, recursive = TRUE)
    log_step(paste("Removed raw GDC download directory:", raw_download_dir))
  }
} else {
  log_step("Cached RNA-seq data found, loading from disk...")
  se <- readRDS(se_path)
}

# --- Clinical data -----------------------------------------------------------
if (!file.exists(clinical_path)) {
  log_step(paste("Querying GDC for", PROJECT, "clinical data..."))
  clinical <- GDCquery_clinic(project = PROJECT, type = "clinical")
  saveRDS(clinical, clinical_path)
} else {
  log_step("Cached clinical data found, loading from disk...")
  clinical <- readRDS(clinical_path)
}

# colnames(clinical)                  # uncomment to explore available variables
# table(clinical$vital_status)

log_step(sprintf("Downloaded: %d samples, %d genes.", ncol(se), nrow(se)))
log_step("01_download_data.R finished.")
