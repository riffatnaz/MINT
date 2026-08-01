# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# run_pipeline.R
# ---------------------------------------------------------------------------
# Driver script -- runs the full pipeline end to end. Adapt to a different
# TCGA project, gene of interest, or thresholds by editing config.R only;
# nothing below needs to change.
#
# Each step also runs on its own once the previous step's cached .rds
# files exist in DATA_DIR, so any stage can be re-run in isolation.
#
# The report step (07) renders report.Rmd from whatever
# is already in RESULTS_DIR and opens it in the default browser. Safe to
# comment out if you don't want that.
# ---------------------------------------------------------------------------

steps <- c(
  "R/01_download_data.R",
  "R/02_preprocess_data.R",
  "R/03_deseq2_analysis.R",
  "R/04_functional_enrichment.R",
  "R/05_gsea_analysis.R",
  "R/06_survival_analysis.R",
  "R/07_generate_report.R"   # automatic HTML report + browser launch
)

for (step in steps) {
  cat("\n==============================\n")
  cat("Running:", step, "\n")
  cat("==============================\n")
  source(step)
}
