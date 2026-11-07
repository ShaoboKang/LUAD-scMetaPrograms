# =============================================================================
# Project Configuration File
# =============================================================================
# Before running any script, please modify BASE_DIR to match your local path.
# All scripts assume the following directory structure under BASE_DIR:
#
# BASE_DIR/
# ├── sc_data/               # Raw scRNA-seq Seurat objects (.rds)
# ├── st_datasets/           # Spatial transcriptomics Seurat objects (.rda)
# ├── bulk/                  # Bulk RNA-seq data
# │   ├── TCGA_GEO7.rdata    # Combined TCGA + 7 GEO datasets
# │   ├── stage_vaild_data/  # Stage validation cohorts
# │   └── sur_vaild_data/    # Survival validation cohorts
# ├── msigdb_v2024.1.Hs_GMTs/# MSigDB gene sets
# └── Result_1 ~ Result_8/   # Output directories (will be created if missing)
#
# =============================================================================

# ---- User-configurable paths ----
BASE_DIR <- path.expand(Sys.getenv("LUAD_ROOT", "~/LUAD_project"))
CODE_DIR <- normalizePath(Sys.getenv("LUAD_CODE_ROOT", getwd()), mustWork = TRUE)  # Path to this code repository

# ---- Derived paths (do not modify below unless you know what you're doing) ----
SC_DATA_DIR      <- file.path(BASE_DIR, "sc_data")
ST_DATA_DIR      <- file.path(BASE_DIR, "st_datasets")
BULK_DIR         <- file.path(BASE_DIR, "bulk")
MSIGDB_DIR       <- file.path(BASE_DIR, "msigdb_v2024.1.Hs_GMTs")

# Sub-directories for intermediate outputs (mirrors original Result_1 ~ Result_8)
RESULT_1_DIR     <- file.path(BASE_DIR, "Result_1")
RESULT_2_DIR     <- file.path(BASE_DIR, "Result_2")
RESULT_3_DIR     <- file.path(BASE_DIR, "Result_3")
RESULT_4_DIR     <- file.path(BASE_DIR, "Result_4")
RESULT_5_DIR     <- file.path(BASE_DIR, "Result_5")
RESULT_6_DIR     <- file.path(BASE_DIR, "Result_6")
RESULT_7_DIR     <- file.path(BASE_DIR, "Result_7")
RESULT_8_DIR     <- file.path(BASE_DIR, "Result_8")

# Helper function to create result subdirectories if they don't exist
init_result_dirs <- function() {
  dirs <- c(RESULT_1_DIR, RESULT_2_DIR, RESULT_3_DIR, RESULT_4_DIR,
            RESULT_5_DIR, RESULT_6_DIR, RESULT_7_DIR, RESULT_8_DIR)
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    rd <- file.path(d, "results_data")
    rf <- file.path(d, "results_figure")
    if (!dir.exists(rd)) dir.create(rd, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(rf)) dir.create(rf, recursive = TRUE, showWarnings = FALSE)
  }
}

# Initialize directories on source
init_result_dirs()

message("[config.R] Loaded. BASE_DIR = ", BASE_DIR)

# English aliases for locally named input folders/files. Override without renaming data.
DIAGNOSTIC_DIR <- path.expand(Sys.getenv("LUAD_DIAGNOSTIC_DIR",
  file.path(BULK_DIR, "\u8bca\u65ad\u6570\u636e")))
CLINICAL_FILE <- path.expand(Sys.getenv("LUAD_CLINICAL_FILE",
  file.path(BULK_DIR, "7\u5957GEO\u548cTCGA\u4e34\u5e8a\u4fe1\u606f\u7528\u4e8e\u591a\u56e0\u7d20fenx.rds")))
MODEL_DIR <- file.path(RESULT_8_DIR, "Mime_TF_selection")
SHAP_DIR <- file.path(RESULT_8_DIR, "SHAP_results")
MP_LEGACY_ORDER <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4",
                     "MP_6", "MP_5", "MP_3", "MP_9")
MP_STATE_NAMES <- c("G2M", "G1S", "NFkB", "MHC-II", "Translation",
                    "OXPHOS", "Adhesion", "EMT", "Glycolysis")
SIGNATURE_TFS <- c("FOXM1", "TFDP1", "FOSL2", "KLF4", "CEBPG", "ATF1", "POU6F1", "XBP1")
FINAL_MODEL <- "StepCox[forward] + Enet[\u03b1=0.7]"
