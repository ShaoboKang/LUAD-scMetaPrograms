# ============================================================================
# Script: 05_immune_celltype_infiltration
# Module: 08_clinical
# Description: Focused CIBERSORT cell-type infiltration comparison for
#   selected immune cell types across datasets.
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

#' Run focused immune cell-type infiltration analysis
#'
#' @param cibersort_path Path to Cibersort RDS from 04_immune_infiltration.
#' @param riskscore_path Path to Model_riskscore RDS.
#' @param output_fig_dir Output directory for figures.
#' @return List of ggplot objects.
run_analysis <- function(
    cibersort_path = file.path(RESULT_8_DIR, "results_data", "Cibersort.rds"),
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "Model_riskscore.rds"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(IOBR)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  library(patchwork)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading CIBERSORT and risk score data...")
  cibersort_list <- readRDS(cibersort_path)
  model_riskscore <- readRDS(riskscore_path)

  plot_list <- list()
  for (i in names(cibersort_list)) {
    message("Processing dataset: ", i)

    cibersort <- cibersort_list[[i]] %>%
      select(c("ID", "T_cells_CD4_memory_activated", "Mast_cells_resting"))

    cibersort$Group <- model_riskscore[[i]]$group[
      match(cibersort$ID, rownames(model_riskscore[[i]]))
    ]
    cibersort$Group <- factor(cibersort$Group, levels = c("Low", "High"))

    res_long <- melt(cibersort, id.vars = c("ID", "Group"),
                     variable.name = "CellType", value.name = "Proportion")

    p <- ggplot(res_long, aes(x = CellType, y = Proportion, fill = Group)) +
      geom_boxplot(outlier.shape = NA, lwd = 0.3, width = 0.7,
                   position = position_dodge(width = 0.75)) +
      ggtitle(i) +
      geom_point(aes(color = Group),
                 position = position_jitterdodge(dodge.width = 0.75,
                                                 jitter.width = 0.15),
                 size = 0.7, alpha = 0.5) +
      stat_compare_means(aes(group = Group), method = "wilcox.test",
                         label = "p.signif", label.y.npc = "top", size = 4) +
      theme_bw(base_size = 12) +
      labs(x = NULL, y = "Score") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 12)
      ) +
      scale_fill_manual(values = c("Low" = "#4575b4", "High" = "#FF7F0E")) +
      scale_color_manual(values = c("Low" = "#4575b4", "High" = "#FF7F0E"))

    plot_list[[i]] <- p
  }

  combined_plot <- wrap_plots(plot_list, ncol = 4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")

  ggsave(
    filename = file.path(output_fig_dir, "06_cibersort_datasets_combined_celltype.pdf"),
    plot = combined_plot,
    width = 15,
    height = 10
  )

  message("Focused immune cell-type analysis complete.")
  return(plot_list)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
