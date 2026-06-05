# ============================================================================
# Script: 10_stage_riskgroup_distribution
# Module: 08_clinical
# Description: Stage distribution by risk group in TCGA with chi-square test.
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

#' Plot stage distribution across risk groups
#'
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param riskscore_path Path to Model_riskscore RDS.
#' @param output_fig_dir Directory for figure output.
#' @return ggplot object.
run_analysis <- function(
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "Model_riskscore.rds"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(dplyr)
  library(ggplot2)
  library(scales)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading bulk data and risk scores...")
  load(bulk_path)
  names(clin_list) <- names(expr_list)

  model_riskscore <- readRDS(file = riskscore_path)

  tcga_clinical <- clin_list$TCGA
  tcga_clinical$RS <- model_riskscore$TCGA$RS
  tcga_clinical$group <- ifelse(tcga_clinical$RS > median(tcga_clinical$RS), "High", "Low")
  tcga_clinical$group <- factor(tcga_clinical$group, levels = c("Low", "High"))
  tcga_clinical <- tcga_clinical[which(tcga_clinical$stage != "'--"), ]

  # Calculate stage frequencies within each risk group
  stage_counts <- tcga_clinical %>%
    group_by(group, stage) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(prop = count / sum(count))

  stage_colors <- c(
    "Stage I"   = "#4DBBD5FF",
    "Stage II"  = "pink",
    "Stage III" = "#F39B7FFF",
    "Stage IV"  = "#E64B35FF"
  )

  p <- ggplot(stage_counts, aes(x = group, y = prop, fill = stage)) +
    geom_col(position = "fill", width = 0.6) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    labs(
      x = "Risk group",
      y = "Proportion",
      title = "Stage distribution by risk group in TCGA"
    ) +
    scale_fill_manual(values = stage_colors) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
    )

  # Chi-square test
  contingency_table <- table(tcga_clinical$group, tcga_clinical$stage)
  chisq_test <- chisq.test(contingency_table)
  p_value <- chisq_test$p.value
  p_label <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3)))

  p <- p + annotate("text", x = 1.5, y = 0.95, label = p_label, size = 5)
  print(p)

  pdf(file = file.path(output_fig_dir, "Stage_distribution_by_risk_group_in_TCGA.pdf"))
  print(p)
  dev.off()

  message("Stage distribution analysis complete.")
  return(p)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
