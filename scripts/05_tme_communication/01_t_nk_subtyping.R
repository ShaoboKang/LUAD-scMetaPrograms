# =============================================================================
# Script: 01_t_nk_subtyping
# Module: 05_tme_communication
# Description: Subtyping of T/NK cells in LUAD scRNA-seq data
# =============================================================================
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggsci)
library(patchwork)
library(reshape2)
library(ggpubr)
library(cowplot)
library(tibble)

#' Load T/NK cell data from the annotated Seurat object
#'
#' @param input_path Path to the annotated Seurat object.
#' @return Seurat object subsetted for T/NK cells.
load_t_nk_data <- function(input_path) {
  readRDS(file = input_path) %>%
    subset(cell_type == "T/NK")
}

#' Process T/NK cells with normalization, PCA, Harmony, and clustering
#'
#' @param obj Seurat object.
#' @return Processed Seurat object with UMAP and clusters.
process_t_nk <- function(obj) {
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()

  # Direct integration UMAP
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

  # Batch correction with Harmony
  set.seed(123)
  obj <- RunHarmony(obj, "Study")
  obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

  # Clustering
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

  return(obj)
}

#' Annotate T/NK cell subtypes based on cluster membership
#'
#' @param obj Seurat object with seurat_clusters.
#' @return Seurat object with T_NK_subtype column.
annotate_t_nk_subtypes <- function(obj) {
  obj$T_NK_subtype <- case_when(
    obj$seurat_clusters %in% c(1, 17, 18) ~ "FCGR3A+ NK",
    obj$seurat_clusters %in% c(12) ~ "FCGR3A- NK",
    obj$seurat_clusters %in% c(5, 6, 7, 9, 10, 13) ~ "CD8+ Cytotoxic",
    obj$seurat_clusters %in% c(3, 15) ~ "CD8+ Exhaustion",
    obj$seurat_clusters %in% c(11) ~ "CD4+ Exhaustion",
    obj$seurat_clusters %in% c(4) ~ "CD4+ Treg",
    obj$seurat_clusters %in% c(0, 2, 8, 14) ~ "CD4+ Naive",
    obj$seurat_clusters %in% c(16) ~ "Epi"
  )
  return(obj)
}

#' Create a dot plot for marker genes
#'
#' @param obj Seurat object.
#' @param features List of marker genes.
#' @param group_by Grouping variable.
#' @return ggplot object.
make_dot_plot <- function(obj, features, group_by) {
  DotPlot(object = obj, features = features, group.by = group_by) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(title = "Expression(%)", order = 3),
      color = guide_colorbar(title = "Average Expression")
    ) +
    scale_color_gradientn(values = seq(0, 1, 0.2), colours = c("#f7f6fb", "#532b83"))
}

#' Plot cell proportion boxplots comparing Normal vs Tumor
#'
#' @param obj Seurat object with subtype and sample metadata.
#' @param subtype_col Column name for subtypes.
#' @param output_path Output PDF path.
#' @param nrow Number of rows for plot grid.
#' @param ncol Number of columns for plot grid.
plot_cell_ratios <- function(obj, subtype_col, output_path, nrow = 2, ncol = 4) {
  groups <- levels(obj@meta.data[[subtype_col]])

  Cellratio <- prop.table(table(obj@meta.data[[subtype_col]], obj$Sample), margin = 2) %>%
    data.frame()
  Cellratio <- dcast(Cellratio, Var2 ~ Var1, value.var = "Freq")
  Cellratio <- column_to_rownames(Cellratio, var = "Var2")

  Cellratio$stage <- as.character(obj$stage[match(rownames(Cellratio), obj$Sample)])
  Cellratio$stage[Cellratio$stage %in% c("Early", "Advance")] <- "Tumor"
  Cellratio$stage <- factor(Cellratio$stage, levels = c("Normal", "Tumor"))

  # Set comparison groups
  group_levels <- levels(factor(Cellratio$stage))
  comp <- combn(group_levels, 2)
  my_comparisons <- list()
  for (j in 1:ncol(comp)) {
    my_comparisons[[j]] <- comp[, j]
  }

  sce_groups <- groups
  col <- c("#1F77B4", "#D62728")

  plist <- list()
  for (group_ in sce_groups) {
    cellper_ <- Cellratio %>% select(one_of(c("stage", group_)))
    colnames(cellper_) <- c("stage", "percent")
    cellper_$percent <- as.numeric(cellper_$percent)

    pb1 <- ggboxplot(
      cellper_,
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
    pb1 <- pb1 + theme(axis.line = element_line(colour = "black")) +
      theme(axis.title.x = element_blank())
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

  pdf(file = output_path, width = 15, height = 8)
  print(plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, align = "vh"))
  dev.off()
}

#' Run T/NK cell subtyping analysis
#'
#' @return Invisible NULL.
run_t_nk_subtyping <- function() {
  input_path <- file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds")
  fig_dir <- file.path(RESULT_5_DIR, "results_figure")
  data_dir <- file.path(RESULT_5_DIR, "results_data")

  # Load and process
  T_NK <- load_t_nk_data(input_path)
  T_NK <- process_t_nk(T_NK)

  # Marker list
  vaildmaker <- list(
    "NK" = c("KLRD1", "GNLY", "KLRB1", "AREG", "FCGR3A"),
    "T" = c("CD3D", "CD3E", "CD3G"),
    "CD8+" = c("CD8A", "CD8B"),
    "CD4+" = c("CD4"),
    "Cytotoxic" = c("GZMA", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7"),
    "Exhaustion" = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4"),
    "Tregs" = c("IL2RA", "FOXP3", "IKZF2"),
    "Naive" = c("TCF7", "SELL", "LEF1", "CCR7"),
    "Epithelial" = c("EPCAM", "MUC1", "KRT19", "KRT18")
  )

  # Plot integration comparison
  p1 <- DimPlot(T_NK, reduction = "umap.unintegrated", group.by = "Study") + theme(legend.position = "none")
  p2 <- DimPlot(T_NK, reduction = "umap.harmony", group.by = c("Study"))
  print(p1 + p2)

  # Annotate subtypes
  T_NK <- annotate_t_nk_subtypes(T_NK)
  T_NK$seurat_clusters <- factor(T_NK$seurat_clusters, levels = rev(c(1, 17, 18, 12, 5, 6, 7, 9, 10, 13, 3, 15, 11, 4, 0, 2, 8, 14, 16)))

  # Cluster marker plot
  pdf(file = file.path(fig_dir, "01_T_cluster_marker.pdf"), width = 10, height = 6)
  print(make_dot_plot(T_NK, vaildmaker, "seurat_clusters"))
  dev.off()

  T_NK$T_NK_subtype <- factor(T_NK$T_NK_subtype, levels = rev(c(
    "FCGR3A+ NK", "FCGR3A- NK", "CD8+ Cytotoxic", "CD8+ Exhaustion",
    "CD4+ Exhaustion", "CD4+ Treg", "CD4+ Naive", "Epi"
  )))

  # Subtype marker plot
  pdf(file = file.path(fig_dir, "02_T_subtype_marker.pdf"), width = 10, height = 6)
  print(make_dot_plot(T_NK, vaildmaker, "T_NK_subtype"))
  dev.off()

  T_NK$T_NK_subtype <- factor(T_NK$T_NK_subtype, levels = c(
    "FCGR3A+ NK", "FCGR3A- NK", "CD8+ Cytotoxic", "CD8+ Exhaustion",
    "CD4+ Exhaustion", "CD4+ Treg", "CD4+ Naive", "Epi"
  ))

  # UMAP plots
  p3 <- DimPlot(T_NK, reduction = "umap.harmony", label = TRUE, group.by = c("seurat_clusters")) + theme(legend.position = "none")
  p4 <- DimPlot(T_NK, reduction = "umap.harmony", label = TRUE, group.by = c("T_NK_subtype")) + scale_color_lancet() + theme(legend.position = "none")

  pdf(file = file.path(fig_dir, "03_T_cluster_subtype_umap.pdf"), width = 13, height = 6)
  print(p3 + p4)
  dev.off()

  saveRDS(T_NK, file = file.path(data_dir, "01_T_NK_anno.rds"))

  # Remove ambiguous clusters
  T_NK <- subset(T_NK, !T_NK_subtype %in% "Epi")
  p5 <- DimPlot(T_NK, reduction = "umap.harmony", label = TRUE, group.by = c("T_NK_subtype")) + scale_color_lancet() + theme(legend.position = "none")

  T_NK$T_NK_subtype <- factor(T_NK$T_NK_subtype, levels = rev(c(
    "FCGR3A+ NK", "FCGR3A- NK", "CD8+ Cytotoxic", "CD8+ Exhaustion",
    "CD4+ Exhaustion", "CD4+ Treg", "CD4+ Naive"
  )))

  vaildmaker <- list(
    "NK" = c("KLRD1", "GNLY", "KLRB1", "AREG", "FCGR3A"),
    "T" = c("CD3D", "CD3E", "CD3G"),
    "CD8+" = c("CD8A", "CD8B"),
    "CD4+" = c("CD4"),
    "Cytotoxic" = c("GZMA", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7"),
    "Exhaustion" = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4"),
    "Tregs" = c("IL2RA", "FOXP3", "IKZF2"),
    "Naive" = c("TCF7", "SELL", "LEF1", "CCR7")
  )

  p6 <- make_dot_plot(T_NK, vaildmaker, "T_NK_subtype")

  pdf(file = file.path(fig_dir, "04_T_cluster_subtype_final.pdf"), width = 15, height = 6)
  print(p5 + p6 + plot_layout(widths = c(1, 2)))
  dev.off()

  saveRDS(T_NK, file = file.path(data_dir, "02_T_NK_anno_final.rds"))

  # Cell ratio boxplots
  T_NK$T_NK_subtype <- factor(T_NK$T_NK_subtype, levels = c(
    "FCGR3A+ NK", "FCGR3A- NK", "CD8+ Cytotoxic", "CD8+ Exhaustion",
    "CD4+ Exhaustion", "CD4+ Treg", "CD4+ Naive"
  ))

  plot_cell_ratios(T_NK, "T_NK_subtype", file.path(fig_dir, "05_T_cluster_subtype_ratio.pdf"))

  message("[INFO] T/NK subtyping completed.")
  invisible(NULL)
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_t_nk_subtyping()
} else {
  message("[INFO] Script loaded. Call run_t_nk_subtyping() to execute.")
}
