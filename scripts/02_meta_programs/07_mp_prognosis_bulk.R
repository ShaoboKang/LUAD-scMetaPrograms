# ============================================================================
# Script: 07_mp_prognosis_bulk
# Module: 02_meta_programs
# Description:
#   Evaluate MP prognostic value across bulk RNA-seq survival cohorts
#   using GSVA/ssGSEA scores and Kaplan-Meier survival curves.
# ============================================================================

# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
  library(GSVA)
  library(tinyarray)
  library(survminer)
  library(survival)

# ---- Helper Functions ----

#' Load MP gene lists.
#' @param rdata_path Path to RData file.
#' @return Named list of MP gene vectors.
load_mp_list <- function(rdata_path) {
  message("Loading MP list from: ", rdata_path)
  load(rdata_path)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_list <- MP_list[old_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))
  mp_list
}

#' Score a survival dataset with ssGSEA.
#' @param file_path Path to RData file.
#' @param mp_list MP gene lists.
#' @return List with gsva_matrix and clinical data frame.
score_survival_dataset <- function(file_path, mp_list) {
  message("Processing: ", basename(file_path))
  load(file_path)
  param <- ssgseaParam(exprData = as.matrix(GSE_expr), geneSets = mp_list)
  gsva_matrix <- gsva(param)
  list(gsva = gsva_matrix, clinical = clinical)
}

#' Generate survival plots for a scored dataset.
#' @param scored_data List from score_survival_dataset.
#' @param mp_list MP gene lists.
#' @param dataset_name Dataset identifier.
#' @return List of ggsurvplot objects.
make_survival_plots <- function(scored_data, mp_list, dataset_name) {
  gsva_matrix <- scored_data$gsva
  clinical <- scored_data$clinical
  surv <- data.frame(
    row.names = colnames(gsva_matrix),
    sample = rownames(clinical),
    event = clinical$vital_status,
    time = clinical$survival_time
  )

  # Divide into high/low groups by median
  for (i in 1:nrow(gsva_matrix)) {
    group <- ifelse(gsva_matrix[i, ] >= quantile(gsva_matrix[i, ])[[3]], "High", "Low")
    surv <- cbind(surv, group)
    colnames(surv)[3 + i] <- rownames(gsva_matrix)[i]
  }

  splots <- list()
  for (o in 1:nrow(gsva_matrix)) {
    fit <- survfit(Surv(time, event) ~ surv[, 3 + o], data = surv)
    splots[[colnames(surv)[3 + o]]] <- ggsurvplot(
      fit, data = surv, palette = c("#E41A1C", "#377EB8"),
      pval = TRUE, size = 1.2, legend = "none", pval.size = 5,
      title = paste(colnames(surv)[3 + o], dataset_name),
      legend.lab = c("High", "Low")
    )
  }
  list(splots = splots, surv_data = surv, gsva = gsva_matrix)
}

# ---- Main Function ----

#' Run MP prognosis bulk analysis.
#' @param rdata_path Path to MP RData file.
#' @param output_fig_dir Directory for output figures.
#' @return List of survival plot lists per dataset.
run_analysis <- function(
    rdata_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_list <- load_mp_list(rdata_path)

  surv_files <- c(
    file.path(BULK_DIR, "sur_vaild_data", "TCGA_sur.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE13213_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE26939_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE30219_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE31210_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE41271_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE72094_sur_data.Rdata"),
    file.path(BULK_DIR, "sur_vaild_data", "GSE11969_sur_data.Rdata")
  )
  dataset_names <- c("TCGA", "GSE13213", "GSE26939", "GSE30219", "GSE31210",
                     "GSE41271", "GSE72094", "GSE11969")

  combined_list <- list()
  for (j in 1:length(surv_files)) {
    scored <- score_survival_dataset(surv_files[j], mp_list)
    res <- make_survival_plots(scored, mp_list, dataset_names[j])
    combined_list <- append(combined_list, res$splots)
  }

  pdf(file = file.path(output_fig_dir, "13_MP_prognosis.pdf"),
      width = 25, height = 20, onefile = FALSE)
  arrange_ggsurvplots(
    combined_list,
    print = TRUE,
    title = "Survival curves",
    ncol = 8, nrow = 9,
    risk.table.height = 0.25
  )
  dev.off()
  message("Saved survival curves to: ", file.path(output_fig_dir, "13_MP_prognosis.pdf"))

  message("MP prognosis analysis complete.")
  invisible(combined_list)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP prognosis analysis: ", e$message)
      stop(e)
    }
  )
}
