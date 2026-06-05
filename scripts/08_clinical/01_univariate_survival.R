# ============================================================================
# Script: 01_univariate_survival
# Module: 08_clinical
# Description: Univariate Cox regression for MP-state-specific TFs across
#   TCGA+7 GEO cohorts, Lasso-Cox feature selection, Kaplan-Meier
#   visualization, and multivariate Cox validation. Generates the 8-TF
#   prognostic signature.
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

#' Univariate Cox regression for a single gene
#'
#' Fits a Cox proportional hazards model for one gene and returns HR, z,
#' p-value, confidence interval, and proportional hazards test p-value.
#'
#' @param gene Gene symbol.
#' @param expr Expression matrix (genes x samples).
#' @param clin Clinical data frame with survival_time and vital_status.
#' @return One-row data frame with Cox results.
gene_cox_result <- function(gene, expr, clin) {
  clin$gene <- t(expr[gene, , drop = FALSE])
  res_cox <- coxph(Surv(survival_time, vital_status) ~ gene, data = clin)
  ph_hypo_multi <- cox.zph(res_cox)$table[1, ]
  cox_summary <- summary(res_cox)
  cox_output <- data.frame(
    gene = gene,
    HR = as.numeric(cox_summary$coefficients[, "exp(coef)"])[1],
    z = as.numeric(cox_summary$coefficients[, "z"])[1],
    pvalue = as.numeric(cox_summary$coefficients[, "Pr(>|z|)"])[1],
    lower = as.numeric(cox_summary$conf.int[, 3][1]),
    upper = as.numeric(cox_summary$conf.int[, 4][1]),
    ph_p = ph_hypo_multi[3]
  )
  return(cox_output)
}

#' Log-rank p-value for a gene set
#'
#' Calculates risk score as mean expression of the gene set and compares
#' high vs low groups by median split using a log-rank test.
#'
#' @param genesets Character vector of gene symbols.
#' @param expr Expression matrix.
#' @param clin Clinical data frame.
#' @return Numeric p-value.
gene_logrank_p <- function(genesets, expr, clin) {
  clin$riskscore <- rowSums(t(expr[genesets, , drop = FALSE]), na.rm = TRUE) / length(genesets)
  breaks <- quantile(clin$riskscore, probs = 1 / 2)
  clin$risk_group <- cut(clin$riskscore, breaks = c(-Inf, breaks, Inf), labels = c("Low", "High"))
  sdiff <- survdiff(Surv(survival_time, vital_status) ~ risk_group, data = clin)
  p_val <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  return(p_val)
}

#' Kaplan-Meier analysis with LASSO coefficients
#'
#' Computes a weighted risk score, dichotomizes by median, and plots a
#' Kaplan-Meier curve with log-rank p-value.
#'
#' @param coefficients Named numeric vector of LASSO coefficients.
#' @param expr Expression matrix.
#' @param clin Clinical data frame.
#' @param time_col Name of survival time column.
#' @param status_col Name of event status column.
#' @param plot_title Plot title string.
#' @return ggsurvplot object.
gene_km_analysis <- function(coefficients, expr, clin,
                             time_col = "survival_time",
                             status_col = "vital_status",
                             plot_title = "Kaplan-Meier Analysis") {
  common_genes <- intersect(names(coefficients), rownames(expr))
  if (length(common_genes) == 0) stop("No matching genes found!")

  clin$riskscore <- colSums(expr[common_genes, , drop = FALSE] * coefficients[common_genes])
  cutoff <- median(clin$riskscore, na.rm = TRUE)
  clin$risk_group <- ifelse(clin$riskscore > cutoff, "High", "Low")
  clin$risk_group <- factor(clin$risk_group, levels = c("Low", "High"))

  sdiff <- survdiff(Surv(survival_time, vital_status) ~ risk_group, data = clin)
  p_val <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)

  km_fit <- survfit(Surv(survival_time, vital_status) ~ risk_group, data = clin)
  km_plot <- ggsurvplot(km_fit, data = clin,
                        title = plot_title,
                        font.x = 14, font.y = 14,
                        font.tickslab = 12,
                        pval = TRUE, pval.size = 4.5,
                        size = 1.2,
                        palette = c("#377EB8", "#E41A1C"),
                        legend = c(0.8, 0.8),
                        legend.labs = c("Low", "High"),
                        legend.title = "Group",
                        tables.theme = theme_cleantable())

  return(km_plot)
}

# ---- Main Function ----

#' Run univariate survival analysis and build prognostic signature
#'
#' @param tf_path Path to top5 TF RDS.
#' @param bulk_path Path to combined TCGA+GEO bulk data.
#' @param output_data_dir Output directory for RDS files.
#' @param output_fig_dir Output directory for figures.
#' @return List with Cox results, LASSO object, and risk score.
run_analysis <- function(
    tf_path = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    output_data_dir = file.path(RESULT_8_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(survival)
  library(survminer)
  library(glmnet)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(forestplot)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading TF names and bulk data...")
  tf_names <- readRDS(tf_path)
  tf_names <- rev(levels(tf_names$Topic))

  tf_all <- gsub("\\(.*\\)", "", tf_names)
  tf_class <- c(rep("MP_1", 5), rep("MP_2", 9), rep("MP_3", 9), rep("MP_4", 9),
                rep("MP_5", 6), rep("MP_7", 7), rep("MP_8", 6), rep("MP_9", 5))
  tf_all_df <- data.frame(Gene = tf_all, CellType = tf_class)

  load(bulk_path)

  # Univariate Cox across all datasets
  message("Running univariate Cox regression...")
  tfs_cox_list <- list()
  for (j in seq_along(expr_list)) {
    tfs_cox <- NULL
    for (i in tf_all) {
      tryCatch({
        tf_cox <- gene_cox_result(i, expr_list[[j]], clin_list[[j]])
        tfs_cox <- rbind(tfs_cox, tf_cox)
      }, error = function(e) {
        message("Error for gene ", i, " in dataset ", names(expr_list)[j])
      })
    }
    tfs_cox$dataset <- names(expr_list)[j]
    tfs_cox_list[[j]] <- tfs_cox
  }
  tfs_cox_all <- do.call(rbind, tfs_cox_list)

  # TCGA significant TFs forest plot
  tfs_cox_sig <- filter(tfs_cox_list[[1]], pvalue < 0.05)
  tfs_cox_sig <- merge(tfs_cox_sig, tf_all_df, by.x = "gene", by.y = "Gene")
  tfs_cox_sig$CellType <- factor(tfs_cox_sig$CellType,
                                 levels = c("MP_1", "MP_2", "MP_3", "MP_4", "MP_7", "MP_8", "MP_9"))
  tfs_cox_sig <- tfs_cox_sig[order(tfs_cox_sig$CellType), ]

  tfs_cox_sig$pvalue_rounded <- sapply(tfs_cox_sig$pvalue, function(x) {
    formatted <- format(round(x, 2), nsmall = 2)
    if (formatted == "0.00") {
      return("0.00")
    }
    return(formatted)
  })

  tabletext <- cbind(tfs_cox_sig$gene, tfs_cox_sig$pvalue_rounded)

  pdf(file.path(output_fig_dir, "01_forestplot.pdf"), width = 5, height = 8, onefile = FALSE)
  forestplot(
    labeltext = tabletext,
    mean = tfs_cox_sig$HR,
    lower = tfs_cox_sig$lower,
    upper = tfs_cox_sig$upper,
    xlog = FALSE,
    main = "Forest Plot of Genes",
    xticks = c(0.5, 1, 1.5, 2),
    col = fpColors(box = "#d43232", line = "#417bab", summary = "royalblue"),
    boxsize = 0.4,
    lwd.ci = 2.6
  )
  dev.off()

  # Z-score boxplot across datasets
  merged_data <- merge(tfs_cox_all, tf_all_df, by.x = "gene", by.y = "Gene")
  p_zscore <- ggplot(merged_data, aes(x = z, y = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
    geom_vline(xintercept = c(-2.58, -1.96, 1.96, 2.58),
               linetype = "dashed", color = "red", size = 0.5) +
    geom_boxplot(show.legend = FALSE, outlier.shape = NA, size = 0.5) +
    geom_jitter(height = 0.15, alpha = 1, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    theme_minimal() +
    theme(
      axis.line.y = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    labs(x = "<--- CoxPH Zscore --->")
  print(p_zscore)

  # Log-rank p-values for individual TFs in TCGA
  message("Calculating log-rank p-values...")
  tfs_lg <- NULL
  for (i in tf_all) {
    tryCatch({
      tf_lg <- gene_logrank_p(i, expr_list[[1]], clin_list[[1]])
      names(tf_lg) <- i
      tfs_lg <- c(tfs_lg, tf_lg)
    }, error = function(e) {
      message("Error for gene ", i)
    })
  }
  message("Significant TFs by log-rank: ")
  print(names(tfs_lg)[tfs_lg < 0.05])

  # Lasso-Cox feature selection (100 iterations)
  message("Running Lasso-Cox selection...")
  set.seed(521)
  genesets <- tfs_cox_sig$gene[-1]
  lasso_result <- list()
  for (j in 1:100) {
    set.seed(j)
    expr <- expr_list[[1]]
    clin <- clin_list[[1]]
    geneset_filter <- intersect(rownames(expr), genesets)
    gene_matrix <- t(as.matrix(expr[geneset_filter, ]))
    response <- Surv(clin$survival_time, clin$vital_status)

    la_eq <- glmnet(gene_matrix, response, family = "cox")
    mod_cv <- cv.glmnet(gene_matrix, response, family = "cox")
    best_model <- glmnet(gene_matrix, response, family = "cox", lambda = mod_cv$lambda.min)
    lasso_result[[j]] <- coef(best_model)
  }
  lasso <- as.matrix(coef(best_model))
  lasso <- lasso[lasso != 0, ]
  genesets <- names(lasso)

  # Kaplan-Meier plots for all datasets
  message("Generating Kaplan-Meier plots...")
  result <- list()
  for (i in seq_along(expr_list)) {
    result[[i]] <- gene_km_analysis(
      coefficients = lasso,
      expr = expr_list[[i]],
      clin = clin_list[[i]],
      plot_title = names(expr_list)[i]
    )
  }

  pdf(file.path(output_fig_dir, "02_model.pdf"), width = 16, height = 8, onefile = FALSE)
  arrange_ggsurvplots(result, ncol = 4, nrow = 2)
  dev.off()

  message("Selected genes in Lasso: ")
  print(tfs_cox_sig[tfs_cox_sig$gene %in% names(lasso), ])

  saveRDS(lasso, file = file.path(output_data_dir, "01_lasso.rds"))

  # ROC / Multivariate Cox using TCGA risk score
  message("Running multivariate Cox validation...")
  common_genes <- intersect(names(lasso), rownames(expr_list[[1]]))
  risk_score <- as.data.frame(colSums(expr_list[[1]][common_genes, , drop = FALSE] * lasso[common_genes]))
  risk_score$sample <- rownames(risk_score)
  colnames(risk_score)[1] <- "riskScore"

  surv <- clin_list[[1]][, c("vital_status", "survival_time", "stage", "gender", "age")]
  surv$sample <- rownames(surv)

  risk_score_cli <- risk_score %>%
    inner_join(surv, by = "sample")
  risk_score_cli$sample <- NULL
  risk_score_cli <- risk_score_cli[risk_score_cli$stage != "'--", ]
  risk_score_cli$stage <- ifelse(risk_score_cli$stage %in% c("Stage I ", "Stage II"), "Early", "Advance")

  multicox <- coxph(Surv(time = survival_time, event = vital_status) ~ ., data = risk_score_cli)
  multisum <- summary(multicox)

  gene_vars <- c("riskScore", "stage", "gender", "age")
  hr <- multisum$coefficients[, 2]
  l95ci <- multisum$conf.int[, 3]
  h95ci <- multisum$conf.int[, 4]
  pvalue <- multisum$coefficients[, 5]
  multiresult <- data.frame(
    gene = gene_vars,
    HR = hr,
    L95CI = l95ci,
    H95CI = h95ci,
    pvalue = pvalue
  )

  a <- multiresult[, c("gene", "pvalue")]
  forestplot(
    a,
    labeltext = a,
    mean = multiresult$HR,
    lower = multiresult$L95CI,
    upper = multiresult$H95CI,
    xlog = FALSE,
    main = "Forest Plot of Genes",
    xticks = c(0.5, 1, 1.5, 2),
    col = fpColors(box = "#d43232", line = "#417bab", summary = "royalblue"),
    boxsize = 0.1,
    lwd.ci = 2.6
  )

  message("Univariate survival analysis complete.")
  return(list(cox_all = tfs_cox_all, lasso = lasso, multiresult = multiresult))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
