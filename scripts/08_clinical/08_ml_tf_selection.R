# ============================================================================
# Script: 08_ml_tf_selection
# Module: 08_clinical
# Description: Machine-learning-based TF prognostic signature selection using
#   Mime1 across batch-corrected TCGA+7 GEO cohorts.
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

# ---- Helper Functions ----

#' Extract gene names from an ML model object
#'
#' Tries multiple slots to retrieve selected features.
#'
#' @param ml_res Single ML model result object.
#' @param model_name Name of the model.
#' @return Character vector of gene names or NULL.
extract_model_genes <- function(ml_res, model_name) {
  gene_model <- ml_res[["xvar.names"]]

  if (is.null(gene_model)) {
    gene_model <- ml_res[["xnames"]]
  }

  if (is.null(gene_model)) {
    coef_min <- coef(ml_res, s = "lambda.min")
    if (!is.numeric(coef_min)) {
      coef_min <- coef_min[which(coef_min[, 1] != 0), ]
      gene_model <- names(coef_min)
    } else {
      coef_min <- coef_min[which(coef_min != 0)]
      gene_model <- names(coef_min)
    }
  }

  if (is.null(gene_model)) {
    gene_model <- ml_res[["fit"]][["var.names"]]
  }

  if (is.null(gene_model)) {
    gene_model <- ml_res[["var.names"]]
  }

  if (is.null(gene_model)) {
    coef <- ml_res[["fit"]][["feature.scores"]]
    if (!is.null(coef)) {
      coef <- as.data.frame(as.matrix(coef))
      coef$gene <- rownames(coef)
      coef <- coef[which(coef[, 1] != 0), ]
      gene_model <- coef$gene
    }
  }

  if (is.null(gene_model)) {
    gene_model <- rownames(ml_res[["pp"]])
  }

  if (is.null(gene_model) && model_name == "Ridge") {
    coef_min <- coef(ml_res[["fit"]],
                     s = ml_res[["cv.fit"]]$lambda.min)
    coef_min <- coef_min[which(coef_min[, 1] != 0), ]
    gene_model <- names(coef_min)
  }

  if (is.null(gene_model) && model_name == "RSF + SuperPC") {
    feature_scores <- ml_res[["RSF + SuperPC"]][[1]][["feature.scores"]]
    feature_scores <- feature_scores[which(feature_scores != 0)]
    gene_model <- names(feature_scores)
  }

  return(gene_model)
}

# ---- Main Function ----

#' Run ML-based TF selection with batch correction
#'
#' @param tf_path Path to top5 TF RDS.
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param ml_func_path Path to modified ML function script.
#' @param output_dir Base output directory for Mime results.
#' @param output_data_dir Directory for result RDS files.
#' @param output_fig_dir Directory for result figures.
#' @return List with ML result, model risk scores, and gene lists.
run_analysis <- function(
    tf_path = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    ml_func_path = file.path(CODE_DIR, "scripts", "08_clinical", "utils",
                             "my_ML.Dev.Prog.Sig_modified_functions.R"),
    output_dir = file.path(RESULT_8_DIR, "Mime_TF_selection"),
    output_data_dir = file.path(RESULT_8_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(Mime1)
  library(doParallel)
  library(tibble)
  library(randomForestSRC)
  library(survminer)
  library(glmnet)
  library(sva)
  library(ggplot2)
  library(dplyr)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading TF names and bulk data...")
  tf_names <- readRDS(tf_path)
  tf_names <- gsub("\\(.*\\)", "", unique(as.character(tf_names$Topic)))

  load(bulk_path)
  names(clin_list) <- names(expr_list)

  common_genes <- Reduce(intersect, lapply(expr_list, rownames))

  # Combine expression matrices for batch correction
  plot_mat <- cbind(
    expr_list$TCGA[common_genes, ], expr_list$GSE11969[common_genes, ],
    expr_list$GSE13213[common_genes, ], expr_list$GSE26939[common_genes, ],
    expr_list$GSE30219[common_genes, ], expr_list$GSE31210[common_genes, ],
    expr_list$GSE41271[common_genes, ], expr_list$GSE72094[common_genes, ]
  )

  dataset_names <- c("TCGA", "GSE11969", "GSE13213", "GSE26939",
                     "GSE30219", "GSE31210", "GSE41271", "GSE72094")
  batch_info <- c(
    rep("TCGA", ncol(expr_list$TCGA)),
    rep("GSE11969", ncol(expr_list$GSE11969)),
    rep("GSE13213", ncol(expr_list$GSE13213)),
    rep("GSE26939", ncol(expr_list$GSE26939)),
    rep("GSE30219", ncol(expr_list$GSE30219)),
    rep("GSE31210", ncol(expr_list$GSE31210)),
    rep("GSE41271", ncol(expr_list$GSE41271)),
    rep("GSE72094", ncol(expr_list$GSE72094))
  )

  message("Running ComBat batch correction...")
  expr_combat <- ComBat(dat = plot_mat, batch = batch_info)

  data_list_combat <- lapply(dataset_names, function(ds) {
    samples <- batch_info == ds
    as.data.frame(t(expr_combat[, samples]))
  })
  names(data_list_combat) <- dataset_names

  data_list <- list()
  for (i in names(data_list_combat)) {
    clin_list[[i]] <- rownames_to_column(clin_list[[i]], var = "ID")
    clin_list[[i]] <- clin_list[[i]][, c("ID", "survival_time", "vital_status")]
    colnames(clin_list[[i]]) <- c("ID", "OS.time", "OS")
    data_list[[i]] <- cbind(clin_list[[i]], data_list_combat[[i]])
  }

  message("Sourcing modified ML function...")
  source(ml_func_path)

  message("Running ML.Dev.Prog.Sig_modified...")
  res <- ML.Dev.Prog.Sig_modified(
    train_data = data_list$TCGA,
    list_train_vali_Data = data_list,
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = 0.05,
    candidate_genes = tf_names,
    mode = "all",
    nodesize = 5,
    seed = 520
  )

  saveRDS(res, file.path(output_dir, "Mime_res_batch_corrected.rds"))

  model_riskscore <- res[["riskscore"]][["StepCox[forward] + Enet[α=0.7]"]]
  saveRDS(model_riskscore, file = file.path(output_data_dir, "Model_riskscore.rds"))

  tcga_riskscore <- res[["riskscore"]][["StepCox[forward] + Enet[α=0.7]"]][["TCGA"]]
  saveRDS(tcga_riskscore, file = file.path(output_data_dir, "TCGA_riskscore.rds"))

  # C-index distribution
  pdf(file.path(output_dir, "cindex_dis_all_batch_corrected.pdf"),
      width = 13, height = 16)
  cindex_dis_all(res, validate_set = names(data_list)[-1],
                 order = names(data_list), width = 0.35)
  dev.off()

  # Survival plots
  survplot <- vector("list", 8)
  for (i in seq_along(data_list)) {
    survplot[[i]] <- rs_sur(
      res,
      model_name = "StepCox[forward] + Enet[α=0.7]",
      dataset = names(data_list)[i],
      median.line = "hv",
      cutoff = 0.5,
      conf.int = TRUE,
      xlab = "Day",
      pval.coord = c(1000, 0.9)
    )
  }

  pdf(file.path(output_dir, "survplot_batch_corrected.pdf"),
      width = 20, height = 10)
  aplot::plot_list(gglist = survplot, ncol = 4)
  dev.off()

  # AUC calculations
  all_auc_1y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = data_list[["TCGA"]],
    inputmatrix.list = data_list,
    mode = "all",
    AUC_time = 1,
    auc_cal_method = "KM"
  )
  all_auc_3y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = data_list[["TCGA"]],
    inputmatrix.list = data_list,
    mode = "all",
    AUC_time = 3,
    auc_cal_method = "KM"
  )
  all_auc_5y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = data_list[["TCGA"]],
    inputmatrix.list = data_list,
    mode = "all",
    AUC_time = 5,
    auc_cal_method = "KM"
  )

  auc_dis_all(all_auc_5y,
              dataset = names(data_list),
              validate_set = names(data_list)[-1],
              order = names(data_list),
              width = 0.35,
              year = 5)

  auc_1y <- roc_vis(
    all_auc_1y,
    model_name = "StepCox[forward] + Enet[α=0.7]",
    dataset = names(data_list),
    order = names(data_list),
    anno_position = c(0.65, 0.55),
    year = 1
  )
  auc_3y <- roc_vis(
    all_auc_3y,
    model_name = "StepCox[forward] + Enet[α=0.7]",
    dataset = names(data_list),
    order = names(data_list),
    anno_position = c(0.65, 0.55),
    year = 3
  )
  auc_5y <- roc_vis(
    all_auc_5y,
    model_name = "StepCox[forward] + Enet[α=0.7]",
    dataset = names(data_list),
    order = names(data_list),
    anno_position = c(0.65, 0.55),
    year = 5
  )

  pdf(file.path(output_dir, "AUC_batch_corrected.pdf"), width = 20, height = 7)
  print(auc_1y + auc_3y + auc_5y)
  dev.off()

  # Extract genes from all models
  message("Extracting genes from all ML models...")
  model_names <- names(res$ml.res)
  genelist <- list()

  for (i in model_names) {
    message("Processing model: ", i)
    gene_model <- extract_model_genes(res[["ml.res"]][[i]], i)
    genelist[[i]] <- gene_model
  }

  # Plot gene inclusion frequency
  gene_df <- stack(genelist)
  colnames(gene_df) <- c("Gene", "Method")

  gene_freq <- gene_df %>%
    group_by(Gene) %>%
    summarise(Frequency = n()) %>%
    arrange(Frequency)
  gene_freq$Gene <- factor(gene_freq$Gene, levels = gene_freq$Gene)

  p_freq <- ggplot(gene_freq, aes(x = Frequency, y = Gene)) +
    geom_segment(aes(x = 0, xend = Frequency, y = Gene, yend = Gene),
                 color = "#C49A6C", linewidth = 0.8) +
    geom_point(aes(size = Frequency),
               color = "#E6A15B", fill = "#E6A15B", shape = 21, alpha = 0.9) +
    scale_size(range = c(3, 10)) +
    theme_classic() +
    labs(x = NULL, y = NULL, size = "Frequency") +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  print(p_freq)

  ggsave(
    filename = file.path(output_fig_dir, "Key_gene_frequency_plot.pdf"),
    plot = p_freq,
    width = 5,
    height = 4
  )

  # Extract coefficients from the best model
  fit_cv <- res$ml.res["StepCox[forward] + Enet[α=0.7]"][[1]]
  coef_min <- coef(fit_cv, s = "lambda.min")
  coef_min <- fit_cv$glmnet.fit$beta[, fit_cv$index["min", "Lambda"]]
  coef_nonzero <- coef_min[coef_min != 0, drop = FALSE]

  message("ML TF selection complete.")
  return(list(res = res, genelist = genelist, coef_nonzero = coef_nonzero))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
