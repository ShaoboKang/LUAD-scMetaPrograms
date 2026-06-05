#' Ligand-Receptor Pair Expression in Cancer Regions
# Script: 09_lr_cancer_region
#'
#' Visualizes spatial transcriptomics data for cancer versus non-cancer regions
#' across selected samples for specific ligand-receptor gene pairs.
#'
#' @return None. PDF figures are written to the result figure directory.
#' @author Script refactored from original analysis pipeline.
#' @export
run_lr_cancer_region <- function() {
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
  .libPaths(c("~/R/SeuratV4", .libPaths()))
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(ggpubr)
  library(lemon)
  library(patchwork)
  library(ggvenn)

  CARD_list <- readRDS(file = file.path(RESULT_6_DIR, "results_data", "CARD_list.rds"))

  # Define gene pairs to analyze
  gene_pairs <- c("MDK", "NCL")
  gene_pairs <- c("SCGB3A2", "MARCO")
  gene_pairs <- c("ANGPTL4", "CDH5")
  gene_pairs <- c("ANGPTL4", "CDH11")
  gene_pairs <- c("PLAU", "PLAUR")

  plots_list <- list()
  for (i in c(
    "seurat_ST8_AIS", "seurat_ST3_MIA", "seurat_ST1_IAC",
    "seurat_ST5_AIS", "seurat_ST6_MIA", "seurat_ST2_IAC"
  )) {
    CARD_list[[i]]$subregion <- ifelse(
      CARD_list[[i]]$subregion == "Cancer", "Cancer", "nonCancer"
    )

    p_subregion <- SpatialPlot(CARD_list[[i]], group.by = "subregion") +
      scale_fill_manual(values = c("#6A3D9A", "#F2F2F2")) +
      theme(legend.position = "right")

    # Fetch gene expression data
    mecom_expr <- FetchData(CARD_list[[i]], vars = gene_pairs, layer = "data")
    mecom_expr$subregion <- CARD_list[[i]]$subregion

    # Create grouping based on gene pair
    gene1 <- gene_pairs[1]
    gene2 <- gene_pairs[2]

    group_col_name <- paste0(gene1, "_", gene2, "_group")

    mecom_expr <- mecom_expr %>%
      mutate(
        !!group_col_name := case_when(
          subregion == "nonCancer" ~ "nonCancer",
          !!sym(gene1) > 0 & !!sym(gene2) > 0 ~ paste0(gene1, "+", gene2, "+"),
          !!sym(gene1) > 0 & !!sym(gene2) == 0 ~ paste0(gene1, "+"),
          !!sym(gene1) == 0 & !!sym(gene2) > 0 ~ paste0(gene2, "+"),
          !!sym(gene1) == 0 & !!sym(gene2) == 0 ~ paste0(gene1, "-", gene2, "-"),
          TRUE ~ "nonCancer"
        )
      )

    # Convert grouping to factor with defined order
    group_levels <- c(
      paste0(gene1, "+", gene2, "+"),
      paste0(gene1, "+"),
      paste0(gene2, "+"),
      paste0(gene1, "-", gene2, "-"),
      "nonCancer"
    )

    mecom_expr[[group_col_name]] <- factor(mecom_expr[[group_col_name]], levels = group_levels)

    # Add to Seurat object
    CARD_list[[i]][[group_col_name]] <- mecom_expr[[group_col_name]]

    # Color scheme for this gene pair
    group_colors <- c("#DC0000FF", "#3C5488FF", "#8491B4FF", "#91D1C2FF", "#F2F2F2")

    # Plot grouping spatial map
    p_group <- SpatialPlot(CARD_list[[i]], group.by = group_col_name) +
      scale_fill_manual(values = group_colors) +
      theme(legend.position = "right")

    plots_list[[i]] <- p_subregion / p_group
  }

  # Create final layout
  final_layout <- wrap_plots(plots_list, ncol = 6) +
    plot_layout(widths = rep(1, 6))

  # Save as PDF
  pdf(
    file.path(RESULT_5_DIR, "results_figure", paste0(gene_pairs[1], "_", gene_pairs[2], "_pairs_ST.pdf")),
    width = 25, height = 6
  )
  print(final_layout)
  dev.off()
}

if (!interactive()) {
  run_lr_cancer_region()
}
