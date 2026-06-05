# =============================================================================
# Script: 02_b_subtyping
# Module: 05_tme_communication
# Description: Subtyping of B cells in LUAD scRNA-seq data
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

#' Load B cell data from the annotated Seurat object
#'
#' @param input_path Path to the annotated Seurat object.
#' @return Seurat object subsetted for B cells.
load_b_data <- function(input_path) {
  readRDS(file = input_path) %>%
    subset(cell_type %in% c("B"))
}

#' Process B cells with normalization, PCA, Harmony, and clustering
#'
#' @param obj Seurat object.
#' @return Processed Seurat object with UMAP and clusters.
process_b_cells <- function(obj) {
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()

  # Direct integration UMAP
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

  # Batch correction with Harmony
  set.seed(123)
  obj <- RunHarmony(obj, "Sample")
  obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

  # Clustering
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

  return(obj)
}

#' Find top 5 DEGs per cluster after removing ambiguous clusters
#'
#' @param obj Seurat object.
#' @param exclude_clusters Clusters to exclude.
#' @return Data frame of top 5 markers per cluster.
find_top5_degs <- function(obj, exclude_clusters) {
  obj_filter <- subset(obj, !seurat_clusters %in% exclude_clusters)
  a <- FindAllMarkers(obj_filter, min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
  a$DEG_pct <- a$pct.1 - a$pct.2
  a_top5 <- a %>%
    group_by(cluster) %>%
    arrange(desc(DEG_pct), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()
  return(a_top5)
}

#' Annotate B cell subtypes based on cluster membership
#'
#' @param obj Seurat object with seurat_clusters.
#' @return Seurat object with B_subtype_subtype column.
annotate_b_subtypes <- function(obj) {
  obj$B_subtype_subtype <- case_when(
    obj$seurat_clusters %in% c(0, 2, 6, 8) ~ "JUND+ B",
    obj$seurat_clusters %in% c(1, 4) ~ "GNB2L1+ B",
    obj$seurat_clusters %in% c(3) ~ "C19ORF43+ B",
    obj$seurat_clusters %in% c(5) ~ "H3+ B",
    obj$seurat_clusters %in% c(9) ~ "FCRL4+ B",
    obj$seurat_clusters %in% c(10) ~ "IFIT3+ B",
    obj$seurat_clusters %in% c(12) ~ "RGS13+ B",
    obj$seurat_clusters %in% c(11) ~ "T",
    obj$seurat_clusters %in% c(7, 13) ~ "Myeloid"
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

#' Run B cell subtyping analysis
#'
#' @return Invisible NULL.
run_b_subtyping <- function() {
  input_path <- file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds")
  fig_dir <- file.path(RESULT_5_DIR, "results_figure")
  data_dir <- file.path(RESULT_5_DIR, "results_data")

  # Load and process
  B_subtype <- load_b_data(input_path)
  B_subtype <- process_b_cells(B_subtype)

  # Plot integration comparison
  p1 <- DimPlot(B_subtype, reduction = "umap.unintegrated", group.by = "Study") + theme(legend.position = "none")
  p2 <- DimPlot(B_subtype, reduction = "umap.harmony", group.by = c("Study"))
  print(p1 + p2)

  # Find top DEGs after removing ambiguous clusters
  a_top5 <- find_top5_degs(B_subtype, c("7", "11", "13"))

  # Marker list
  vaildmaker <- list(
    "B" = c("CD79A", "MS4A1"),
    "Top5_DEGs_in_B_cluster" = c(unique(a_top5$gene)),
    "T" = c("CD3D", "CD3E", "CD3G", "CD8A"),
    "Myeloid" = c("CD68", "LYZ", "AIF1", "CTSD")
  )

  # Annotate subtypes
  B_subtype <- annotate_b_subtypes(B_subtype)
  B_subtype$seurat_clusters <- factor(B_subtype$seurat_clusters, levels = rev(c(0, 2, 6, 8, 1, 4, 3, 5, 9, 10, 12, 11, 7, 13)))

  # Cluster marker plot
  pdf(file = file.path(fig_dir, "06_B_cluster_marker.pdf"), width = 10, height = 6)
  print(make_dot_plot(B_subtype, vaildmaker, "seurat_clusters"))
  dev.off()

  B_subtype$B_subtype_subtype <- factor(B_subtype$B_subtype_subtype, levels = rev(c(
    "JUND+ B", "GNB2L1+ B", "C19ORF43+ B", "H3+ B", "FCRL4+ B",
    "IFIT3+ B", "RGS13+ B", "T", "Myeloid"
  )))

  # Subtype marker plot
  pdf(file = file.path(fig_dir, "07_B_subtype_marker.pdf"), width = 10, height = 6)
  print(make_dot_plot(B_subtype, vaildmaker, "B_subtype_subtype"))
  dev.off()

  B_subtype$B_subtype_subtype <- factor(B_subtype$B_subtype_subtype, levels = c(
    "JUND+ B", "GNB2L1+ B", "C19ORF43+ B", "H3+ B",
    "FCRL4+ B", "IFIT3+ B", "RGS13+ B", "T", "Myeloid"
  ))

  custom_colors <- c(
    "#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF",
    "#6A6599FF", "#80796BFF", "#adb6b699", "#adb6b699"
  )

  # UMAP plots
  p3 <- DimPlot(B_subtype, reduction = "umap.harmony", label = TRUE, group.by = c("seurat_clusters")) + theme(legend.position = "none")
  p4 <- DimPlot(B_subtype, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = c("B_subtype_subtype")) +
    scale_color_manual(values = custom_colors) + theme(legend.position = "none")

  pdf(file = file.path(fig_dir, "08_B_cluster_subtype_umap.pdf"), width = 13, height = 6)
  print(p3 + p4)
  dev.off()

  saveRDS(B_subtype, file = file.path(data_dir, "03_B_subtype_anno.rds"))

  # Remove ambiguous clusters
  B_subtype <- subset(B_subtype, !B_subtype_subtype %in% c("T", "Myeloid"))
  p5 <- DimPlot(B_subtype, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = c("B_subtype_subtype")) +
    scale_color_manual(values = custom_colors) + theme(legend.position = "none")

  B_subtype$B_subtype_subtype <- factor(B_subtype$B_subtype_subtype, levels = rev(c(
    "JUND+ B", "GNB2L1+ B", "C19ORF43+ B", "H3+ B",
    "FCRL4+ B", "IFIT3+ B", "RGS13+ B"
  )))

  vaildmaker <- list(
    "B" = c("CD79A", "MS4A1"),
    "Top5_DEGs_in_B_cluster" = c(unique(a_top5$gene))
  )

  p6 <- make_dot_plot(B_subtype, vaildmaker, "B_subtype_subtype")

  pdf(file = file.path(fig_dir, "09_B_cluster_subtype_final.pdf"), width = 15, height = 6)
  print(p5 + p6 + plot_layout(widths = c(1, 2)))
  dev.off()

  saveRDS(B_subtype, file = file.path(data_dir, "04_B_subtype_anno_final.rds"))

  # Cell ratio boxplots
  B_subtype$B_subtype_subtype <- factor(B_subtype$B_subtype_subtype, levels = c(
    "JUND+ B", "GNB2L1+ B", "C19ORF43+ B", "H3+ B",
    "FCRL4+ B", "IFIT3+ B", "RGS13+ B"
  ))

  plot_cell_ratios(B_subtype, "B_subtype_subtype", file.path(fig_dir, "10_B_cluster_subtype_ratio.pdf"))

  message("[INFO] B cell subtyping completed.")
  invisible(NULL)
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_b_subtyping()
} else {
  message("[INFO] Script loaded. Call run_b_subtyping() to execute.")
}
