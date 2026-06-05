# ============================================================================
# Script: 06_drug_sensitivity
# Module: 08_clinical
# Description: oncoPredict drug sensitivity prediction using CTRP2 training
#   data across TCGA+7 GEO cohorts.
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

# ---- Main Function ----

#' Run drug sensitivity prediction with oncoPredict
#'
#' @param ctrp2_expr_path Path to CTRP2 expression RDS.
#' @param ctrp2_res_path Path to CTRP2 response RDS.
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param output_data_dir Directory for predictions.
#' @return List of prediction matrices per dataset.
run_analysis <- function(
    ctrp2_expr_path = file.path(RESULT_8_DIR, "DataFiles", "Training Data",
                                "CTRP2_Expr (TPM, not log transformed).rds"),
    ctrp2_res_path = file.path(RESULT_8_DIR, "DataFiles", "Training Data", "CTRP2_Res.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    output_data_dir = file.path(RESULT_8_DIR, "results_data")
) {
  library(oncoPredict)
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)

  message("Loading CTRP2 training data...")
  ctrp2_expr <- readRDS(ctrp2_expr_path)
  ctrp2_expr <- log2(ctrp2_expr + 1)

  ctrp2_res <- readRDS(ctrp2_res_path)

  message("Loading bulk data...")
  load(bulk_path)

  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)

  pred_ic50_results <- list()

  for (i in names(expr_list)) {
    message("Predicting drug sensitivity for ", i, "...")
    ds_dir <- file.path(output_data_dir, paste0("ic50_", i))
    if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

    old_wd <- getwd()
    setwd(ds_dir)
    on.exit(setwd(old_wd), add = TRUE)

    pred_ic50_results[[i]] <- calcPhenotype(
      trainingExprData = ctrp2_expr,
      trainingPtype = ctrp2_res,
      testExprData = as.matrix(expr_list[[i]]),
      batchCorrect = "eb",
      minNumSamples = 20,
      printOutput = TRUE,
      removeLowVaryingGenes = 0.2,
      removeLowVaringGenesFrom = "homogenizeData"
    )

    setwd(old_wd)
  }

  saveRDS(pred_ic50_results,
          file = file.path(output_data_dir, "pred_ic50_results.rds"))

  message("Drug sensitivity prediction complete.")
  return(pred_ic50_results)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
