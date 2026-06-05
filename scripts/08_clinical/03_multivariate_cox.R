# ============================================================================
# Script: 03_multivariate_cox
# Module: 08_clinical
# Description: Multivariate Cox regression forest plots for TCGA using
#   stage or T/N/M stratification.
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

#' Build a forest plot from a multivariate Cox result data frame
#'
#' @param df Data frame with Variable, HR, Lower95, Upper95, pValue columns.
#' @param x_lab X-axis label.
#' @return Invisible NULL; plot is printed.
plot_multivariate_forest <- function(df, x_lab = "Hazard Ratios") {
  p_vals <- df$pValue
  p_text <- ifelse(
    p_vals >= 1e-4,
    sprintf("%.4f", p_vals),
    sprintf("%.2e", p_vals)
  )

  tabletext <- cbind(
    c(NA, as.character(df$Variable)),
    c("HR", sprintf("%.2f", df$HR)),
    c("95% CI", paste0(sprintf("%.2f", df$Lower95), "–", sprintf("%.2f", df$Upper95))),
    c("p.value", p_text)
  )

  hl <- list(
    "1" = gpar(lwd = 2, col = "black"),
    "2" = gpar(lwd = 1.5, col = "black")
  )
  hl[[as.character(nrow(df) + 2)]] <- gpar(lwd = 2, col = "black")

  xt <- round(c(min(df$Lower95) - 0.1, 1, max(df$Upper95) + 0.1), 1)

  forestplot(
    labeltext = tabletext,
    mean = c(NA, df$HR),
    lower = c(NA, df$Lower95),
    upper = c(NA, df$Upper95),
    graph.pos = 4,
    graphwidth = unit(0.40, "npc"),
    fn.ci_norm = "fpDrawCircleCI",
    col = fpColors(
      box = "#00A087FF",
      lines = "black",
      zero = "#D0D1E6"
    ),
    boxsize = 0.35,
    lwd.ci = 1.5,
    ci.vertices = TRUE,
    ci.vertices.height = 0.1,
    zero = 1,
    lwd.zero = 1.5,
    xticks = xt,
    lwd.xaxis = 2,
    xlab = x_lab,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.3),
      ticks = gpar(cex = 1.1),
      xlab = gpar(cex = 1.3),
      title = gpar(cex = 1.5)
    ),
    hrzl_lines = hl,
    lineheight = unit(0.9, "cm"),
    colgap = unit(0.35, "cm"),
    mar = unit(rep(1.2, 4), "cm"),
    new_page = FALSE
  )
}

#' Summarize a coxph fit into a tidy data frame
#'
#' @param fit A coxph model fit.
#' @param var_labels Character vector of display names matching coef rows.
#' @return Tidy data frame.
summarize_cox_fit <- function(fit, var_labels) {
  su <- summary(fit)
  ci_tab <- as.data.frame(su$conf.int[, c("exp(coef)", "lower .95", "upper .95")])
  p_tab <- su$coefficients[, "Pr(>|z|)", drop = FALSE]
  df <- cbind(
    Variable = rownames(ci_tab),
    HR = ci_tab$`exp(coef)`,
    Lower95 = ci_tab$`lower .95`,
    Upper95 = ci_tab$`upper .95`,
    pValue = p_tab[, 1]
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$HR <- as.numeric(df$HR)
  df$Lower95 <- as.numeric(df$Lower95)
  df$Upper95 <- as.numeric(df$Upper95)
  df$pValue <- as.numeric(df$pValue)
  df$Variable <- var_labels
  df$Variable <- factor(df$Variable, levels = var_labels)
  df <- df[order(df$Variable), ]
  return(df)
}

# ---- Main Function ----

#' Run multivariate Cox regression and generate forest plots
#'
#' @param riskscore_path Path to TCGA risk score RDS.
#' @param clinical_path Path to TCGA clinical RData.
#' @param output_fig_dir Directory for PDF outputs.
#' @return List with both model summaries.
run_analysis <- function(
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "TCGA_riskscore.rds"),
    clinical_path = file.path(BULK_DIR, "TCGA", "TCGA_sur_data.Rdata"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(survival)
  library(survminer)
  library(forestplot)
  library(grid)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading risk score and clinical data...")
  tcga_riskscore <- readRDS(riskscore_path)
  load(clinical_path)
  rm(mRNA_exp)
  gc()

  com_sample <- intersect(rownames(tcga_riskscore), clinical$case_submitter_id)
  tcga_riskscore <- tcga_riskscore[com_sample, , drop = FALSE]
  clinical <- clinical[match(com_sample, clinical$case_submitter_id), ]
  clinical$riskscore <- tcga_riskscore$RS

  clinical <- clinical[, c("case_submitter_id", "days_to_death", "event",
                            "riskscore", "age_at_index", "gender", "stage",
                            "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")]
  colnames(clinical) <- c("Sample_ID", "OS_time", "OS", "Riskscore",
                          "Age", "Gender", "stage", "T_stage", "N_stage", "M_stage")
  clinical$Age <- as.numeric(clinical$Age)
  clinical$stage <- factor(clinical$stage, levels = c("early", "advance"))

  # Multivariate Cox with stage
  message("Fitting multivariate Cox model with stage...")
  fmla <- Surv(OS_time, OS) ~ Riskscore + Age + Gender + stage
  fit <- coxph(fmla, data = clinical)
  df_stage <- summarize_cox_fit(fit, c("Riskscore", "Age", "Gender", "Stage"))

  pdf(file.path(output_fig_dir, "TCGA_multivariate_cox_forestplot.pdf"),
      width = 8, height = 4.5)
  plot_multivariate_forest(df_stage, "Hazard Ratios")
  dev.off()

  # Multivariate Cox with T/N/M
  message("Fitting multivariate Cox model with T/N/M...")
  clinical$T_stage <- gsub("[abc]", "", clinical$T_stage)
  clinical$M_stage <- ifelse(grepl("M1", clinical$M_stage), "M1", "M0")
  clinical$T_stage[clinical$T_stage == "TX"] <- NA
  clinical$N_stage[clinical$N_stage == "NX"] <- NA
  clinical$M_stage[clinical$M_stage == "MX"] <- NA
  clinical <- na.omit(clinical)

  clinical$T_stage <- ifelse(clinical$T_stage %in% c("T1", "T2"), "T1-2", "T3-4")
  clinical$N_stage <- ifelse(clinical$N_stage == "N0", "N0", "N1-3")

  fmla_tnm <- Surv(OS_time, OS) ~ Riskscore + Age + Gender + T_stage + N_stage + M_stage
  fit_tnm <- coxph(fmla_tnm, data = clinical)
  df_tnm <- summarize_cox_fit(fit_tnm, c("Riskscore", "Age", "Gender",
                                         "T_stage", "N_stage", "M_stage"))

  pdf(file.path(output_fig_dir, "TCGA_multivariate_cox_forestplot_TNM.pdf"),
      width = 8, height = 4.5)
  plot_multivariate_forest(df_tnm, "Hazard Ratios")
  dev.off()

  message("Multivariate Cox analysis complete.")
  return(list(stage = df_stage, tnm = df_tnm))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
