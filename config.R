# MINT 1.0.0 -- Modular Integrated Next-gen Transcriptomics 

# Transcriptomic & Differential Expression Analysis Pipeline

# config.R
# ---------------------------------------------------------------------------
# Central configuration. Change values here to adapt the pipeline to a
# different TCGA project, gene of interest, or thresholds -- no other
# script needs to be touched.
# ---------------------------------------------------------------------------

PROJECT     <- "TCGA-BRCA"   # any valid TCGA project code, e.g. "TCGA-LUAD"
DATA_DIR    <- "data"        # cached raw/intermediate data (created if missing)
RESULTS_DIR <- "results"     # all tables and plots are written here

# Sample types to download/keep.
# NOTE: not valid for every project -- TCGA-LAML (a blood cancer) uses
# "Primary Blood Derived Cancer - Peripheral Blood" / "Blood Derived Normal"
# instead. Check available sample types before running on a non solid-tumor
# project.
SAMPLE_TYPES <- c("Primary Tumor", "Solid Tissue Normal")

# Minimum Normal-condition samples before 02_preprocess_data.R warns that
# the Tumor-vs-Normal comparison may be underpowered.
MIN_NORMAL_SAMPLES <- 10

# Differential expression thresholds. Applied throughout: DESeq2 DEG
# calling, functional enrichment (GO/KEGG/Reactome), GSEA, and the report.
PVAL_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# LFC shrinkage method for the MA plot only (deseq2_results.csv/GSEA stay
# unshrunken). "ashr" (recommended, works directly with a contrast),
# "apeglm" (needs a coef name, not a contrast), "normal", or "none" to
# disable shrinkage.
LFC_SHRINK_METHOD <- "ashr"

# Gene used for the illustrative expression boxplot and survival split.
GENE_OF_INTEREST <- "ESR1"

# Functional enrichment reference databases.
ORGANISM_DB       <- "org.Hs.eg.db"
KEGG_ORGANISM     <- "hsa"
REACTOME_ORGANISM <- "human"

# GSEA gene-set collection, using msigdbr's collection/subcollection naming
# (run msigdbr::msigdbr_collections() to see all valid combinations).
# Examples: collection = "H" (hallmark, no subcollection needed);
# collection = "C2", subcollection = "CP:KEGG_MEDICUS" (or "CP:REACTOME",
# "CP:WIKIPATHWAYS"); collection = "C5", subcollection = "GO:BP".
GSEA_COLLECTION    <- "C2"
GSEA_SUBCOLLECTION <- "CP:KEGG_MEDICUS"

# Clinical fields used for survival analysis (TCGA clinical field names).
SURV_TIME_VAR  <- "days_to_last_follow_up"
SURV_EVENT_VAR <- "vital_status"

# Optional clinical covariates for the multivariate Cox model, used only if
# present in the clinical data (missing columns are dropped automatically).
COX_COVARIATES <- c("age_at_index", "ajcc_pathologic_stage")

# Optional: collapse sparse factor levels (e.g. AJCC sub-stages) before Cox
# fitting, to avoid separation from categories with very few
# patients/events. Applied in 02_preprocess_data.R. A covariate with no
# entry here is left as-is; leave as list() to skip, or replace with a
# project's own column/levels to adapt.
COVARIATE_COLLAPSE_MAP <- list(
  ajcc_pathologic_stage = c(
    "Stage IA" = "Stage I", "Stage IB" = "Stage I", "Stage I" = "Stage I",
    "Stage IIA" = "Stage II", "Stage IIB" = "Stage II", "Stage II" = "Stage II",
    "Stage IIIA" = "Stage III", "Stage IIIB" = "Stage III", "Stage IIIC" = "Stage III", "Stage III" = "Stage III",
    "Stage IV" = "Stage IV", "Stage X" = NA
  )
)
