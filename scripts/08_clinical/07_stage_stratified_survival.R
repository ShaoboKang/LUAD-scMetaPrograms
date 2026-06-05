# ============================================================================
# Script: 07_stage_stratified_survival
# Module: 08_clinical
# Description: Stage-stratified Kaplan-Meier survival analysis for TCGA
#   colored by risk group.
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

#' Plot Kaplan-Meier curve for a single stage
#'
#' @param dat Clinical data frame with survival_time, event, and RiskGroup.
#' @param stage_name Stage label to filter.
#' @return A ggplot object.
plot_stage_km <- function(dat, stage_name) {
  sub <- dat %>%
    filter(stage == stage_name) %>%
    filter(!is.na(survival_time), !is.na(event), !is.na(RiskGroup))

  n_low  <- sum(sub$RiskGroup == "Low")
  n_high <- sum(sub$RiskGroup == "High")

  if (n_low < 5 || n_high < 5) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(stage_name, "\nInsufficient samples"),
               size = 6) +
      theme_void()
    return(p)
  }

  fit <- survfit(Surv(survival_time, event) ~ RiskGroup, data = sub)
  cox_fit <- coxph(Surv(survival_time, event) ~ RiskGroup, data = sub)
  cox_sum <- summary(cox_fit)
  hr <- round(cox_sum$coefficients[1, "exp(coef)"], 2)
  p_cox <- cox_sum$coefficients[1, "Pr(>|z|)"]

  p_txt <- ifelse(p_cox < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_cox)))
  hr_txt <- paste0("HR = ", hr)

  legend_labs <- c(
    paste0("Low(n=", n_low, ")"),
    paste0("High(n=", n_high, ")")
  )

  g <- ggsurvplot(
    fit,
    data = sub,
    conf.int = TRUE,
    censor = FALSE,
    risk.table = FALSE,
    surv.median.line = "hv",
    pval = FALSE,
    legend.title = NULL,
    legend.labs = legend_labs,
    palette = c("#9E9E9E", "#B45A5A"),
    ggtheme = theme_classic(base_size = 14),
    title = stage_name
  )

  p <- g$plot +
    annotate(
      "text",
      x = max(sub$survival_time, na.rm = TRUE) * 0.05,
      y = 0.18,
      label = paste(p_txt, hr_txt, sep = "\n"),
      hjust = 0,
      size = 5
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Survival probability") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = c(0.80, 0.86),
      legend.background = element_blank(),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "black", size = 11),
      axis.line = element_line(linewidth = 0.9)
    )

  return(p)
}

# ---- Main Function ----

#' Run stage-stratified survival analysis
#'
#' @param riskscore_path Path to TCGA risk score RDS.
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param output_fig_dir Directory for figure output.
#' @return Combined ggplot object.
run_analysis <- function(
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "TCGA_riskscore.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(ggplot2)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading risk score and clinical data...")
  tcga_riskscore <- readRDS(riskscore_path)
  load(bulk_path)
  names(clin_list) <- names(expr_list)

  tcga_clinical <- clin_list$TCGA
  tcga_clinical$RS <- tcga_riskscore$RS
  tcga_clinical <- tcga_clinical[which(tcga_clinical$stage != "'--"), ]
  tcga_clinical$stage <- factor(tcga_clinical$stage,
                                levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

  cutoff <- median(tcga_clinical$RS, na.rm = TRUE)
  tcga_clinical$RiskGroup <- ifelse(tcga_clinical$RS > cutoff, "High", "Low")
  tcga_clinical$RiskGroup <- factor(tcga_clinical$RiskGroup, levels = c("Low", "High"))

  message("Generating stage-stratified KM plots...")
  p1 <- plot_stage_km(tcga_clinical, "Stage I")
  p2 <- plot_stage_km(tcga_clinical, "Stage II")
  p3 <- plot_stage_km(tcga_clinical, "Stage III")
  p4 <- plot_stage_km(tcga_clinical, "Stage IV")

  p_final <- ggarrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    font.label = list(size = 18, face = "bold")
  )

  ggsave(
    filename = file.path(output_fig_dir, "TCGA_stage_stratified_KM.pdf"),
    plot = p_final,
    width = 10,
    height = 8
  )

  message("Stage-stratified survival analysis complete.")
  return(p_final)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
