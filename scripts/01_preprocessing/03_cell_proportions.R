# ============================================================================
# Script: 03_cell_proportions
# Module: 01_preprocessing
# Description:
#   Calculate and visualize cell type proportions across samples and stages
#   using boxplots and stacked barplots.
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

#' Load annotated Seurat object.
#' @param path Path to .rds file.
#' @return Seurat object.
load_annotated_data <- function(path) {
  message("Loading annotated data from: ", path)
  readRDS(path)
}

#' Compute cell type proportions per sample.
#' @param seurat_obj Seurat object with cell_type and Sample metadata.
#' @return Data frame of proportions.
compute_cell_ratios <- function(seurat_obj) {
  cellratio <- prop.table(table(seurat_obj$cell_type, seurat_obj$Sample), margin = 2) %>%
    data.frame()
  cellratio <- dcast(cellratio, Var2 ~ Var1, value.var = "Freq")
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

#' Plot boxplots of cell proportions by stage.
#' @param cellratio Data frame with stage and proportion columns.
#' @param sce_groups Character vector of cell type names to plot.
#' @param my_comparisons List of comparison pairs.
#' @param col Color palette.
#' @return List of ggplot objects.
plot_proportion_boxplots <- function(cellratio, sce_groups, my_comparisons, col) {
  plist <- list()
  for (group_ in sce_groups) {
    cellper <- cellratio %>% select(one_of(c("stage", group_)))
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

#' Plot stacked barplot of cell type proportions.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_proportion_barplot <- function(seurat_obj, output_path, width = 14, height = 7) {
  count_table <- table(seurat_obj$cell_type, seurat_obj$Sample) %>% as.data.frame()
  colnames(count_table) <- c("cell_type", "sample", "count")
  count_table <- count_table %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

  df_wide <- count_table %>%
    pivot_wider(
      id_cols = sample,
      names_from = cell_type,
      values_from = proportion
    )
  df_wide$stage <- seurat_obj$stage[match(df_wide$sample, seurat_obj$Sample)]
  colnames(df_wide)[3] <- "TNK"
  df_wide <- df_wide %>%
    arrange(
      stage,
      Epithelial, TNK, B, Plasma, Myeloid, Mast, Fibroblast, Endothelial
    )

  sample_order <- df_wide %>% pull(sample)
  count_table <- count_table %>%
    mutate(sample = factor(sample, levels = sample_order))

  custom_colors <- c(
    "#C6307C", "#D0AFC4", "#4991C1", "#89558D",
    "#AFC2D9", "#435B95", "#79B99D", "#F5A623"
  )

  p <- ggplot(count_table, aes(x = sample, y = proportion, fill = factor(cell_type))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = custom_colors, name = "Cluster") +
    labs(x = "Sample ID", y = "Proportion", title = "Proportion of cell_type by Sample") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      legend.position = "top"
    )

  pdf(file = output_path, width = width, height = height)
  print(p)
  dev.off()
  message("Saved barplot to: ", output_path)
}

# ---- Main Function ----

#' Run cell proportion analysis.
#' @param input_path Path to annotated Seurat object.
#' @param output_fig_dir Directory for output figures.
#' @return List of plot objects (invisibly).
run_analysis <- function(
    input_path = file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds"),
    output_fig_dir = file.path(RESULT_1_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) {
    dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
  }

  luad_merge <- load_annotated_data(input_path)

  # Compute ratios
  cellratio <- compute_cell_ratios(luad_merge)
  my_comparisons <- make_comparisons(cellratio$stage)

  sce_groups <- levels(luad_merge$cell_type)
  col <- c("#79B99D", "#1F77B4", "#D62728")

  plist <- plot_proportion_boxplots(cellratio, sce_groups, my_comparisons, col)

  pdf(file = file.path(output_fig_dir, "07_celltype_ratio.pdf"), width = 15, height = 8)
  print(plot_grid(plotlist = plist, nrow = 2, ncol = 4, align = "vh"))
  dev.off()
  message("Saved boxplots to: ", file.path(output_fig_dir, "07_celltype_ratio.pdf"))

  plot_proportion_barplot(luad_merge, file.path(output_fig_dir, "08_celltype_ratio_barplot.pdf"))

  invisible(plist)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in cell proportions: ", e$message)
      stop(e)
    }
  )
}
