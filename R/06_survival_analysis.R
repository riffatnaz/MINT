# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# R/06_survival_analysis.R
# ---------------------------------------------------------------------------
# Survival analysis, patients stratified by expression of a gene of
# interest (config.R -> GENE_OF_INTEREST). Swap that one value to analyze
# any other gene without touching this script.
#
# Two complementary approaches:
#   1. Kaplan-Meier -- non-parametric, High/Low expression groups (median
#      split), with a log-rank test.
#   2. Cox proportional hazards -- expression treated as continuous, both
#      univariate and (if clinical covariates are available) multivariate.
# ---------------------------------------------------------------------------

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(rlang)
})

# --- 0. Setup: load data, build analysis frame -------------------------------
setup_ok <- TRUE

# Survival time uses SURV_TIME_VAR (config.R) and falls back to
# days_to_death where that's missing -- TCGA typically records
# days_to_last_follow_up for living patients and days_to_death for
# deceased ones, so a single time column needs both. If SURV_TIME_VAR is
# changed to something other than days_to_last_follow_up, confirm this
# fallback still makes sense for that field.
setup_result <- tryCatch({
  vsd       <- readRDS(file.path(DATA_DIR, "vsd.rds"))
  gene_info <- readRDS(file.path(DATA_DIR, "gene_info.rds"))
  clinical  <- readRDS(file.path(DATA_DIR, "clinical_processed.rds"))
  col_data  <- readRDS(file.path(DATA_DIR, "col_data.rds"))

  gene_id <- gene_info$gene_id[gene_info$gene_name == GENE_OF_INTEREST][1]
  if (is.na(gene_id)) stop("GENE_OF_INTEREST not found in gene_info: ", GENE_OF_INTEREST)

  expr <- SummarizedExperiment::assay(vsd)[gene_id, ]

  df <- data.frame(barcode = names(expr), expression = expr, row.names = NULL) %>%
    filter(col_data[barcode, "condition"] == "Tumor") %>%
    mutate(submitter_id = substr(barcode, 1, 12)) %>%
    inner_join(clinical, by = "submitter_id") %>%
    mutate(
      time  = coalesce(!!sym(SURV_TIME_VAR), days_to_death),
      event = ifelse(!!sym(SURV_EVENT_VAR) == "Dead", 1, 0),
      group = ifelse(expression >= median(expression), "High", "Low")
    ) %>%
    filter(!is.na(time), time >= 0)

  ensure_dir(RESULTS_DIR)
  df
}, error = function(e) {
  log_step(paste("Setup failed --", conditionMessage(e), "-- skipping survival analysis."))
  NULL
})

if (is.null(setup_result)) {
  setup_ok <- FALSE
} else {
  surv_df <- setup_result
}

# --- 1. Kaplan-Meier ----------------------------------------------------------
if (setup_ok) {
  tryCatch({
    fit <- survfit(Surv(time, event) ~ group, data = surv_df)

    km_plot <- ggsurvplot(
      fit, data = surv_df, pval = TRUE, risk.table = TRUE,
      legend.title = GENE_OF_INTEREST,
      title = paste(PROJECT, "-", GENE_OF_INTEREST, "expression and survival")
    )

    ggplot2::ggsave(
      file.path(RESULTS_DIR, "survival_km_plot.png"),
      plot = km_plot$plot, width = 8, height = 7, dpi = 300, bg = "white"
    )

    log_step(sprintf(
      "Kaplan-Meier complete: %d patients (%d High / %d Low).",
      nrow(surv_df), sum(surv_df$group == "High"), sum(surv_df$group == "Low")
    ))
  }, error = function(e) {
    log_step(paste("Kaplan-Meier failed --", conditionMessage(e), "-- continuing to Cox."))
  })
}

# --- 2. Cox proportional hazards ----------------------------------------------

# Strips the covariate name prefix off factor-level terms (e.g.
# "ajcc_pathologic_stageStage III" -> "Stage III").
clean_term_names <- function(terms, covariate_names) {
  for (cov in covariate_names) terms <- sub(paste0("^", cov), "", terms)
  trimws(terms)
}

extract_cox_table <- function(cox_fit, covariate_names = character(0)) {
  coefs     <- summary(cox_fit)$coefficients
  ci        <- summary(cox_fit)$conf.int
  raw_terms <- rownames(coefs)

  # Traces each term back to its source covariate for a separate column.
  matched_cov <- sapply(raw_terms, function(t) {
    hits <- covariate_names[startsWith(t, covariate_names)]
    if (length(hits) == 0) return(NA_character_)
    hits[which.max(nchar(hits))]
  })

  data.frame(
    covariate = matched_cov,
    term      = clean_term_names(raw_terms, covariate_names),
    HR        = coefs[, "exp(coef)"],
    HR_lower  = ci[, "lower .95"],
    HR_upper  = ci[, "upper .95"],
    p_value   = coefs[, grepl("^Pr", colnames(coefs))],
    row.names = NULL
  )
}

if (setup_ok) {
  tryCatch({
    cox_uni <- coxph(Surv(time, event) ~ expression, data = surv_df)
    cox_uni_table <- extract_cox_table(cox_uni, covariate_names = "expression") %>%
      mutate(model = "Univariate (expression only)")
    write.csv(cox_uni_table, file.path(RESULTS_DIR, "cox_univariate.csv"), row.names = FALSE)
    log_step("Univariate Cox regression complete.")
  }, error = function(e) {
    log_step(paste("Univariate Cox failed --", conditionMessage(e), "-- continuing to multivariate."))
  })

  tryCatch({
    # Only covariates present with >1 distinct value are used, so this
    # adapts across projects without editing the script.
    available_covariates <- COX_COVARIATES[
      COX_COVARIATES %in% colnames(surv_df) &
        sapply(COX_COVARIATES, function(v) dplyr::n_distinct(surv_df[[v]], na.rm = TRUE) > 1)
    ]

    if (length(available_covariates) == 0) {
      log_step("Multivariate Cox skipped (no usable covariates found).")
    } else {
      cox_formula <- reformulate(c("expression", available_covariates), response = "Surv(time, event)")
      cox_multi <- coxph(cox_formula, data = surv_df)

      cox_multi_table <- extract_cox_table(cox_multi, covariate_names = c("expression", available_covariates)) %>%
        mutate(model = "Multivariate (adjusted)")
      write.csv(cox_multi_table, file.path(RESULTS_DIR, "cox_multivariate.csv"), row.names = FALSE)

      # Guard: non-finite HR/CI (separation) would crash ggforest()'s log
      # axis. Should be rare now that sparse levels are pre-collapsed in
      # 02_preprocess_data.R, but kept as a safety net.
      finite_hr <- all(is.finite(summary(cox_multi)$conf.int))

      if (finite_hr) {
        forest_plot <- ggforest(cox_multi, data = surv_df)
        ggplot2::ggsave(
          file.path(RESULTS_DIR, "cox_forest_plot.png"),
          plot = forest_plot, width = 8, height = 5 + 0.4 * length(available_covariates),
          dpi = 300, bg = "white"
        )
        log_step(sprintf(
          "Multivariate Cox complete (adjusted for: %s).",
          paste(available_covariates, collapse = ", ")
        ))
      } else {
        log_step("Multivariate Cox has non-finite HR/CI -- forest plot skipped, see cox_multivariate.csv.")
      }
    }
  }, error = function(e) {
    log_step(paste("Multivariate Cox failed --", conditionMessage(e), "-- univariate results still saved."))
  })
}

log_step("06_survival_analysis.R finished.")
