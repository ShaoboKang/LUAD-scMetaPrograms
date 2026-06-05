# =============================================================================
# Script: 03_myeloid_subtyping
# Module: 05_tme_communication
# Description: Subtyping of myeloid cells in LUAD scRNA-seq data
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

#' Load myeloid cell data from the annotated Seurat object
#'
#' @param input_path Path to the annotated Seurat object.
#' @return Seurat object subsetted for myeloid cells.
load_myeloid_data <- function(input_path) {
  readRDS(file = input_path) %>%
    subset(cell_type %in% c("Myeloid"))
}

#' Process myeloid cells with normalization, PCA, Harmony, and clustering
#'
#' @param obj Seurat object.
#' @return Processed Seurat object with UMAP and clusters.
process_myeloid_cells <- function(obj) {
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

#' Find top 5 DEGs per cluster after removing ambiguous clusters
#'
#' @param obj Seurat object.
#' @param exclude_clusters Clusters to exclude.
#' @return Data frame of top 5 markers per cluster.
find_top5_degs <- function(obj, exclude_clusters) {
  obj_filter <- subset(obj, !seurat_clusters %in% exclude_clusters)
  a <- FindAllMarkers(obj_filter, min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
  a$DEG_pct <- a$pct.1 - a$pct.2
  a$cluster <- factor(a$cluster, levels = c(0, 16, 17, 1, 8, 12, 2, 4, 5, 10, 13, 3, 9, 6, 14, 15, 18, 7))
  a_top5 <- a %>%
    group_by(cluster) %>%
    arrange(desc(DEG_pct), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()
  return(a_top5)
}

#' Annotate myeloid cell subtypes based on cluster membership
#'
#' @param obj Seurat object with seurat_clusters.
#' @return Seurat object with Myeloid_subtype_subtype column.
annotate_myeloid_subtypes <- function(obj) {
  obj$Myeloid_subtype_subtype <- case_when(
    obj$seurat_clusters %in% c(0, 16, 17) ~ "TCEB2+ Mac",
    obj$seurat_clusters %in% c(1, 8, 12) ~ "MCEMP1+ Mac",
    obj$seurat_clusters %in% c(2) ~ "SELENOP+ Mac",
    obj$seurat_clusters %in% c(4) ~ "VEGFA+ Mac",
    obj$seurat_clusters %in% c(5) ~ "SPP1+ Mac",
    obj$seurat_clusters %in% c(10) ~ "CXCL10+ Mac",
    obj$seurat_clusters %in% c(13) ~ "C3+ Mac",
    obj$seurat_clusters %in% c(3) ~ "S100A12+ Mono",
    obj$seurat_clusters %in% c(9) ~ "LILRB2+ Mono",
    obj$seurat_clusters %in% c(6) ~ "CD1C+ DC",
    obj$seurat_clusters %in% c(14) ~ "CLEC9A+ DC",
    obj$seurat_clusters %in% c(15) ~ "LAMP3+ DC",
    obj$seurat_clusters %in% c(18) ~ "Neutrophil",
    obj$seurat_clusters %in% c(7) ~ "Low quality Mac",
    obj$seurat_clusters %in% c(11) ~ "T",
    obj$seurat_clusters %in% c(19) ~ "Plasma"
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
#' @param groups Vector of group names to plot. If NULL, uses all levels.
#' @param nrow Number of rows for plot grid.
#' @param ncol Number of columns for plot grid.
plot_cell_ratios <- function(obj, subtype_col, output_path, groups = NULL, nrow = 2, ncol = 6) {
  if (is.null(groups)) {
    groups <- levels(obj@meta.data[[subtype_col]])
  }

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

  pdf(file = output_path, width = 20, height = 8)
  print(plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, align = "vh"))
  dev.off()
}

#' Run myeloid cell subtyping analysis
#'
#' @return Invisible NULL.
run_myeloid_subtyping <- function() {
  input_path <- file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds")
  fig_dir <- file.path(RESULT_5_DIR, "results_figure")
  data_dir <- file.path(RESULT_5_DIR, "results_data")

  # Load and process
  Myeloid_subtype <- load_myeloid_data(input_path)
  Myeloid_subtype <- process_myeloid_cells(Myeloid_subtype)

  # Plot integration comparison
  p1 <- DimPlot(Myeloid_subtype, reduction = "umap.unintegrated", group.by = "Study") + theme(legend.position = "none")
  p2 <- DimPlot(Myeloid_subtype, reduction = "umap.harmony", group.by = c("Study"))
  print(p1 + p2)

  # Find top DEGs after removing ambiguous clusters
  a_top5 <- find_top5_degs(Myeloid_subtype, c("11", "19"))

  # Marker list
  vaildmaker <- list(
    Myeloid = c("CD68", "LYZ"),
    "Top5_DEGs_in_Myeloid_cluster" = c(unique(a_top5$gene)),
    "T" = c("CD3D", "CD3E", "CD3G"),
    "Plasma" = c("MZB1", "JCHAIN", "IGHG1")
  )

  # Annotate subtypes
  Myeloid_subtype <- annotate_myeloid_subtypes(Myeloid_subtype)
  Myeloid_subtype$seurat_clusters <- factor(
    Myeloid_subtype$seurat_clusters,
    levels = rev(c(0, 16, 17, 1, 8, 12, 2, 4, 5, 10, 13, 3, 9, 6, 14, 15, 18, 7, 11, 19))
  )

  # Cluster marker plot
  pdf(file = file.path(fig_dir, "11_Myeloid_cluster_marker.pdf"), width = 12, height = 6)
  print(make_dot_plot(Myeloid_subtype, vaildmaker, "seurat_clusters"))
  dev.off()

  Myeloid_subtype$Myeloid_subtype_subtype <- factor(
    Myeloid_subtype$Myeloid_subtype_subtype,
    levels = rev(c(
      "TCEB2+ Mac", "MCEMP1+ Mac", "SELENOP+ Mac", "VEGFA+ Mac", "SPP1+ Mac",
      "CXCL10+ Mac", "C3+ Mac", "S100A12+ Mono", "LILRB2+ Mono", "CD1C+ DC",
      "CLEC9A+ DC", "LAMP3+ DC", "Neutrophil", "Low quality Mac", "T", "Plasma"
    ))
  )

  # Subtype marker plot
  pdf(file = file.path(fig_dir, "12_Myeloid_subtype_marker.pdf"), width = 12, height = 6)
  print(make_dot_plot(Myeloid_subtype, vaildmaker, "Myeloid_subtype_subtype"))
  dev.off()

  Myeloid_subtype$Myeloid_subtype_subtype <- factor(
    Myeloid_subtype$Myeloid_subtype_subtype,
    levels = rev(levels(Myeloid_subtype$Myeloid_subtype_subtype))
  )

  custom_colors <- c(
    "#FF000099", "#FF990099", "#FFCC0099", "#00FF0099", "#6699FF99", "#CC33FF99",
    "#99991E99", "#99999999", "#FF00CC99", "#CC000099", "#FFCCCC99", "#FFFF0099",
    "#CCFF0099", "grey", "grey", "grey"
  )

  # UMAP plots
  p3 <- DimPlot(Myeloid_subtype, reduction = "umap.harmony", label = TRUE, group.by = c("seurat_clusters")) + theme(legend.position = "none")
  p4 <- DimPlot(Myeloid_subtype, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = c("Myeloid_subtype_subtype")) +
    scale_color_manual(values = custom_colors) + theme(legend.position = "none")

  pdf(file = file.path(fig_dir, "13_Myeloid_cluster_subtype_umap.pdf"), width = 13, height = 6)
  print(p3 + p4)
  dev.off()

  saveRDS(Myeloid_subtype, file = file.path(data_dir, "05_Myeloid_subtype_anno.rds"))

  # Remove ambiguous clusters
  Myeloid_subtype <- subset(Myeloid_subtype, !Myeloid_subtype_subtype %in% c("Low quality Mac", "T", "Plasma"))
  p5 <- DimPlot(Myeloid_subtype, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = c("Myeloid_subtype_subtype")) +
    scale_color_manual(values = custom_colors) + theme(legend.position = "none")

  Myeloid_subtype$Myeloid_subtype_subtype <- factor(
    Myeloid_subtype$Myeloid_subtype_subtype,
    levels = rev(c(
      "TCEB2+ Mac", "MCEMP1+ Mac", "SELENOP+ Mac", "VEGFA+ Mac", "SPP1+ Mac",
      "CXCL10+ Mac", "C3+ Mac", "S100A12+ Mono", "LILRB2+ Mono", "CD1C+ DC",
      "CLEC9A+ DC", "LAMP3+ DC", "Neutrophil"
    ))
  )

  vaildmaker <- list(
    Myeloid = c("CD68", "LYZ"),
    "Top5_DEGs_in_Myeloid_cluster" = c(unique(a_top5$gene))
  )

  p6 <- make_dot_plot(Myeloid_subtype, vaildmaker, "Myeloid_subtype_subtype")

  pdf(file = file.path(fig_dir, "14_Myeloid_cluster_subtype_final.pdf"), width = 18, height = 6)
  print(p5 + p6 + plot_layout(widths = c(1, 3)))
  dev.off()

  saveRDS(Myeloid_subtype, file = file.path(data_dir, "06_Myeloid_subtype_anno_final.rds"))

  # Cell ratio boxplots (first 12 groups only)
  Myeloid_subtype$Myeloid_subtype_subtype <- factor(
    Myeloid_subtype$Myeloid_subtype_subtype,
    levels = rev(levels(Myeloid_subtype$Myeloid_subtype_subtype))
  )

  ratio_groups <- levels(Myeloid_subtype$Myeloid_subtype_subtype)[1:12]
  plot_cell_ratios(
    Myeloid_subtype, "Myeloid_subtype_subtype",
    file.path(fig_dir, "15_Myeloid_cluster_subtype_ratio.pdf"),
    groups = ratio_groups
  )

  message("[INFO] Myeloid subtyping completed.")
  invisible(NULL)
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_myeloid_subtyping()
} else {
  message("[INFO] Script loaded. Call run_myeloid_subtyping() to execute.")
}
