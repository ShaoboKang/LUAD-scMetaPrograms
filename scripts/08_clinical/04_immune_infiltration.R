# ============================================================================
# Script: 04_immune_infiltration
# Module: 08_clinical
# Description: ESTIMATE and CIBERSORT immune infiltration analysis across
#   TCGA+7 GEO cohorts stratified by risk group.
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

#' Run immune deconvolution and generate boxplots for one method
#'
#' @param expr_list List of expression matrices.
#' @param model_riskscore List of risk score data frames.
#' @param method Deconvolution method ("estimate" or "cibersort").
#' @param output_fig Path to output PDF.
#' @param output_rds Path to output RDS.
#' @param width PDF width.
#' @param height PDF height.
#' @return List of deconvolution results.
run_deconvolution <- function(expr_list, model_riskscore, method,
                               output_fig, output_rds,
                               width = 24, height = 10) {
  library(IOBR)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  library(patchwork)

  res_list <- list()
  plot_list <- list()

  for (i in names(expr_list)) {
    message("Running ", method, " for ", i, "...")

    if (method == "estimate") {
      res <- deconvo_tme(eset = expr_list[[i]], method = "estimate")
      colnames(res) <- gsub("_estimate$", "", colnames(res))
      res <- res %>% select(-TumorPurity)
      y_lab <- "Score"
    } else if (method == "cibersort") {
      res <- deconvo_tme(eset = expr_list[[i]], method = "cibersort",
                        arrays = TRUE, perm = 100)
      res <- res %>%
        select(-c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT"))
      colnames(res) <- gsub("_CIBERSORT$", "", colnames(res))
      y_lab <- "Proportion"
    } else {
      stop("Unsupported method: ", method)
    }

    res$Group <- model_riskscore[[i]]$group[match(res$ID, rownames(model_riskscore[[i]]))]
    res$Group <- factor(res$Group, levels = c("Low", "High"))
    res_list[[i]] <- res

    res_long <- melt(res, id.vars = c("ID", "Group"),
                     variable.name = "CellType", value.name = y_lab)

    p <- ggplot(res_long, aes(x = CellType, y = .data[[y_lab]], fill = Group)) +
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
      labs(x = NULL, y = y_lab) +
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

  ggsave(filename = output_fig, plot = combined_plot,
         width = width, height = height)
  saveRDS(res_list, file = output_rds)

  message(method, " analysis complete.")
  return(res_list)
}

# ---- Main Function ----

#' Run immune infiltration analysis
#'
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param riskscore_path Path to Model_riskscore RDS.
#' @param output_data_dir Output directory for RDS.
#' @param output_fig_dir Output directory for figures.
#' @return List with ESTIMATE and CIBERSORT results.
run_analysis <- function(
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "Model_riskscore.rds"),
    output_data_dir = file.path(RESULT_8_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(dplyr)

  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading bulk data and risk scores...")
  load(bulk_path)
  model_riskscore <- readRDS(riskscore_path)
  model_riskscore <- lapply(model_riskscore, function(df) {
    med <- median(df$RS, na.rm = TRUE)
    df$group <- ifelse(df$RS > med, "High", "Low")
    return(df)
  })

  estimate_list <- run_deconvolution(
    expr_list, model_riskscore, method = "estimate",
    output_fig = file.path(output_fig_dir, "05_Estimate_datasets_combined.pdf"),
    output_rds = file.path(output_data_dir, "Estimate.rds"),
    width = 24, height = 10
  )

  cibersort_list <- run_deconvolution(
    expr_list, model_riskscore, method = "cibersort",
    output_fig = file.path(output_fig_dir, "06_cibersort_datasets_combined.pdf"),
    output_rds = file.path(output_data_dir, "Cibersort.rds"),
    width = 35, height = 15
  )

  message("Immune infiltration analysis complete.")
  return(list(estimate = estimate_list, cibersort = cibersort_list))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
