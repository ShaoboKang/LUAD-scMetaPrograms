# ============================================================================
# Script: 09_mp_ratio_analysis
# Module: 02_meta_programs
# Description:
#   Calculate and visualize MP proportions across samples and stages
#   using boxplots for early vs advance and split stage groups.
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
  library(reshape2)
  library(dplyr)
  library(ggpubr)
  library(cowplot)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(Seurat)
  library(tidyverse)

# ---- Helper Functions ----

#' Load assigned MP Seurat object and standardize MP names.
#' @param rds_path Path to RDS file.
#' @return Seurat object with standardized MP column.
load_mp_data <- function(rds_path) {
  message("Loading MP assigned data from: ", rds_path)
  luad_epi <- readRDS(rds_path)
  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))
  luad_epi
}

#' Compute MP proportions per sample.
#' @param seurat_obj Seurat object with MP and Sample metadata.
#' @return Data frame of proportions.
compute_mp_ratios <- function(seurat_obj) {
  cellratio <- prop.table(table(seurat_obj$MP, seurat_obj$Sample), margin = 2) %>%
    data.frame()
  cellratio <- reshape2::dcast(cellratio, Var2 ~ Var1, value.var = "Freq")
  cellratio <- column_to_rownames(cellratio, var = "Var2")
  cellratio$stage <- seurat_obj$stage[match(rownames(cellratio), seurat_obj$Sample)]
  cellratio
}

#' Generate pairwise comparisons for all stage levels.
#' @param stage_vec Character vector of stage values.
#' @return List of comparison pairs.
make_comparisons <- function(stage_vec) {
  group <- levels(factor(stage_vec))
  comp <- combn(group, 2)
  my_comparisons <- list()
  for (j in 1:ncol(comp)) {
    my_comparisons[[j]] <- comp[, j]
  }
  my_comparisons
}

#' Plot boxplots of MP proportions by stage.
#' @param cellratio Data frame with stage and proportion columns.
#' @param sce_groups Character vector of MP names.
#' @param my_comparisons List of comparison pairs.
#' @param col Color palette.
#' @return List of ggplot objects.
plot_mp_proportion_boxplots <- function(cellratio, sce_groups, my_comparisons, col) {
  plist <- list()
  for (group_ in sce_groups) {
    cellper <- cellratio %>% dplyr::select(one_of(c("stage", group_)))
    colnames(cellper) <- c("stage", "percent")
    cellper$percent <- as.numeric(cellper$percent)

    pb1 <- ggboxplot(
      cellper,
      x = "stage",
      y = "percent",
      color = "stage",
      fill = NULL,
      add = "jitter",
      bxp.errorbar.width = 0.8,
      width = 0.5,
      size = 0.5,
      font.label = list(size = 15),
      palette = col
    ) +
      theme(panel.background = element_blank())
    pb1 <- pb1 + theme(axis.line = element_line(colour = "black"))
    pb1 <- pb1 + theme(axis.title.x = element_blank())
    pb1 <- pb1 + theme(axis.title.y = element_blank()) +
      theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1))
    pb1 <- pb1 + theme(axis.text.y = element_text(size = 15)) +
      ggtitle(group_) +
      theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
    pb1 <- pb1 + theme(legend.position = "NA")
    pb1 <- pb1 + stat_compare_means(
      method = "t.test", hide.ns = FALSE,
      comparisons = my_comparisons, label = "p.adj.format",
      p.adjust.method = "BH", vjust = 0.02, bracket.size = 0.6
    )
    plist[[group_]] <- pb1
  }
  plist
}

#' Save a grid of boxplots to PDF.
#' @param plist List of ggplot objects.
#' @param output_path Output PDF path.
#' @param nrow Number of rows.
#' @param ncol Number of columns.
#' @param width PDF width.
#' @param height PDF height.
save_boxplot_grid <- function(plist, output_path, nrow = 3, ncol = 3,
                              width = 12, height = 12) {
  pdf(file = output_path, width = width, height = height)
  print(plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, align = "vh"))
  dev.off()
  message("Saved boxplot grid to: ", output_path)
}

# ---- Main Function ----

#' Run MP ratio analysis.
#' @param rds_path Path to assigned MP Seurat object.
#' @param output_fig_dir Directory for output figures.
#' @return List of plot lists.
run_analysis <- function(
    rds_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  luad_epi <- load_mp_data(rds_path)

  # Early vs Advance
  cellratio <- compute_mp_ratios(luad_epi)
  my_comparisons <- make_comparisons(cellratio$stage)
  sce_groups <- levels(luad_epi$MP)
  col <- c("#1F77B4", "#D62728")

  plist <- plot_mp_proportion_boxplots(cellratio, sce_groups, my_comparisons, col)
  save_boxplot_grid(plist, file.path(output_fig_dir, "04_MP_tumor_stage_boxplot.pdf"))

  # Split stage (I/II/III/IV)
  cellratio$stage <- luad_epi$AJCC_stage[match(rownames(cellratio), luad_epi$Sample)]
  cellratio$stage <- case_when(
    cellratio$stage %in% c("MIA", "AIS", "IA", "IA3", "IAC", "IB", "IB-IIA") ~ "I",
    cellratio$stage %in% c("IIA", "IIB") ~ "II",
    cellratio$stage %in% c("IIIA", "IIIB", "IIIB ") ~ "III",
    cellratio$stage %in% c("IV") ~ "IV"
  )
  cellratio$stage <- factor(cellratio$stage, levels = c("I", "II", "III", "IV"))

  save_boxplot_grid(plist, file.path(output_fig_dir, "04_MP_tumor_stage_boxplot_split.pdf"))

  message("MP ratio analysis complete.")
  invisible(plist)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP ratio analysis: ", e$message)
      stop(e)
    }
  )
}
