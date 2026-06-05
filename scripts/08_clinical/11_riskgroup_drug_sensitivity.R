# ============================================================================
# Script: 11_riskgroup_drug_sensitivity
# Module: 08_clinical
# Description: Compare predicted drug IC50 between risk groups and generate
#   boxplots for selected drugs.
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

#' Run Wilcoxon test per drug between risk groups
#'
#' @param sen Data frame with drug columns, RS, and group.
#' @param drug_cols Character vector of drug column names.
#' @return Data frame with test statistics.
test_drugs_by_group <- function(sen, drug_cols) {
  test_results <- lapply(drug_cols, function(drug) {
    low_vals <- sen[sen$group == "Low", drug]
    high_vals <- sen[sen$group == "High", drug]
    low_vals <- na.omit(low_vals)
    high_vals <- na.omit(high_vals)

    if (length(low_vals) < 2 || length(high_vals) < 2) {
      return(data.frame(
        drug = drug,
        p_value = NA,
        mean_low = NA,
        mean_high = NA,
        statistic = NA,
        n_low = length(low_vals),
        n_high = length(high_vals)
      ))
    }

    test <- wilcox.test(high_vals, low_vals)
    data.frame(
      drug = drug,
      p_value = test$p.value,
      mean_low = mean(low_vals),
      mean_high = mean(high_vals),
      statistic = test$statistic,
      n_low = length(low_vals),
      n_high = length(high_vals)
    )
  })

  results_df <- bind_rows(test_results)
  results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
  return(results_df)
}

# ---- Main Function ----

#' Run risk-group drug sensitivity comparison
#'
#' @param drug_pred_path Path to DrugPredictions CSV.
#' @param riskscore_path Path to Model_riskscore RDS.
#' @param output_fig_dir Directory for figures.
#' @return List with test results and boxplot.
run_analysis <- function(
    drug_pred_path = file.path(RESULT_8_DIR, "results_data",
                               "calcPhenotype_Output", "DrugPredictions.csv"),
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "Model_riskscore.rds"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading drug predictions and risk scores...")
  sen <- read.csv(drug_pred_path, row.names = 1, header = TRUE, check.names = FALSE)
  model_riskscore <- readRDS(file = riskscore_path)

  sen$RS <- model_riskscore$TCGA$RS
  sen$group <- ifelse(sen$RS > median(sen$RS), "High", "Low")
  sen$group <- factor(sen$group, levels = c("Low", "High"))

  drug_cols <- setdiff(names(sen), c("RS", "group"))

  message("Running group comparisons per drug...")
  results_df <- test_drugs_by_group(sen, drug_cols)

  sig_drugs <- results_df %>% filter(p_value < 0.05) %>% arrange(p_value)
  sig_drugs_fdr <- results_df %>% filter(p_adj < 0.05) %>% arrange(p_adj)

  message("Significant drugs (raw p < 0.05):")
  print(sig_drugs)
  message("Significant drugs (FDR < 0.05):")
  print(sig_drugs_fdr)

  # Boxplots for selected drugs
  drug_long <- sen %>%
    tibble::rownames_to_column("Sample_ID") %>%
    pivot_longer(
      cols = all_of(drug_cols),
      names_to = "Drug",
      values_to = "Predicted_IC50"
    )

  drug_order <- c("nutlin-3", "nilotinib", "bexarotene",
                  "MK-2206", "MK-1775", "paclitaxel")

  plot_df <- drug_long %>%
    mutate(Drug = factor(Drug, levels = drug_order))
  plot_df <- plot_df[!is.na(plot_df$Drug), ]

  pdf(file = file.path(output_fig_dir, "Sensitivity_sen_in_TCGA.pdf"),
      height = 3.5, width = 15)
  p <- ggplot(plot_df, aes(x = group, y = Predicted_IC50, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0) +
    facet_wrap(~ Drug, scales = "free_y", ncol = 6) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = c("Low" = "#4575b4", "High" = "#FF7F0E")) +
    theme_classic() +
    labs(y = "Predicted IC50") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  print(p)
  dev.off()

  message("Risk-group drug sensitivity analysis complete.")
  return(list(test_results = results_df, plot = p))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
