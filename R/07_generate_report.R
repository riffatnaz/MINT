# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/07_generate_report.R
# ---------------------------------------------------------------------------
# Renders report.Rmd into a single self-contained HTML report, embedding
# every table and plot produced by the pipeline plus narrative text, and
# opens it in the default browser when done.
#
# Optional and separate from run_pipeline.R -- run it once earlier stages
# have produced their outputs in RESULTS_DIR:
#
#   source("R/07_generate_report.R")
#
# or from the command line:
#
#   Rscript -e 'source("R/07_generate_report.R")'
#
# All values shown in the report (project, gene, cutoffs, results) come
# from config.R and RESULTS_DIR at the time this script runs -- there is
# no separate copy to keep in sync. To refresh the report after changing
# config.R or re-running an analysis stage, just re-run this script.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to generate the report. Install it with install.packages('rmarkdown').")
}

ensure_dir(RESULTS_DIR)

output_file <- paste0("analysis_report_", gsub("[^A-Za-z0-9_-]", "_", PROJECT), ".html")
output_path <- file.path(RESULTS_DIR, output_file)

# Explicit params list passed to rmarkdown::render() below overrides
# report.Rmd's own YAML `params:` defaults -- those defaults are only used
# if someone hits RStudio's "Knit" button directly on report.Rmd, bypassing
# this script. Coercing types here (character/numeric) protects against a
# stale or hand-edited YAML default of the wrong type ever leaking into a
# numeric comparison inside the report.
report_params <- list(
  project     = as.character(PROJECT),
  results_dir = as.character(RESULTS_DIR),
  gene        = as.character(GENE_OF_INTEREST),
  pval_cutoff = as.numeric(PVAL_CUTOFF),
  lfc_cutoff  = as.numeric(LFC_CUTOFF)
)

log_step(paste(
  "Rendering report to", output_path,
  sprintf("(project=%s, gene=%s, pval_cutoff=%s, lfc_cutoff=%s)",
          report_params$project, report_params$gene,
          report_params$pval_cutoff, report_params$lfc_cutoff)
))

rmarkdown::render(
  input       = "report.Rmd",
  output_file = output_file,
  output_dir  = RESULTS_DIR,
  params      = report_params,
  envir       = new.env(),
  quiet       = TRUE
)

log_step("Report generation complete.")

launch_file(output_path)
log_step(paste("Opened report in default browser:", output_path))
