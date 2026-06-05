# ============================================================================
# Script: 02_cell_annotation
# Module: 01_preprocessing
# Description:
#   Batch-correct clustering with Harmony, annotate major cell types,
#   remove doublet clusters, and save final annotated Seurat object.
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
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)

# ---- Helper Functions ----

#' Load the filtered merged Seurat object.
#' @param path Path to the .rds file.
#' @return Seurat object.
load_filtered_data <- function(path) {
  message("Loading filtered data from: ", path)
  readRDS(path)
}

#' Run standard Seurat workflow: Normalize, FindVariableFeatures, Scale, PCA.
#' @param seurat_obj Seurat object.
#' @param nfeatures Number of variable features. Default 2000.
#' @return Processed Seurat object.
run_standard_workflow <- function(seurat_obj, nfeatures = 2000) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj
}

#' Integrate batches with Harmony.
#' @param seurat_obj Seurat object with PCA computed.
#' @param group_by Column to integrate by. Default "Sample".
#' @return Harmony-integrated Seurat object.
run_harmony_integration <- function(seurat_obj, group_by = "Sample") {
  seurat_obj <- RunHarmony(seurat_obj, group.by = group_by)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
  seurat_obj
}

#' Plot UMAP before and after Harmony integration.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_harmony_comparison <- function(seurat_obj, output_path) {
  p1 <- DimPlot(seurat_obj, reduction = "umap.unintegrated", group.by = "Study") +
    theme(legend.position = "none")
  p2 <- DimPlot(seurat_obj, reduction = "umap.harmony", group.by = "Study")
  pdf(file = output_path, width = 15, height = 6)
  print(p1 + p2)
  dev.off()
  message("Saved harmony comparison to: ", output_path)
}

#' Annotate major cell types from cluster assignments.
#' @param seurat_obj Seurat object with seurat_clusters.
#' @return Seurat object with cell_type column.
annotate_cell_types <- function(seurat_obj) {
  seurat_obj$cell_type <- dplyr::case_when(
    seurat_obj$seurat_clusters %in% c(2, 16, 17, 19, 20, 23) ~ "Epithelial",
    seurat_obj$seurat_clusters %in% c(0, 1, 6, 12, 14, 18, 21, 26, 31) ~ "T/NK",
    seurat_obj$seurat_clusters %in% c(4, 24, 25) ~ "B",
    seurat_obj$seurat_clusters %in% c(13) ~ "Plasma",
    seurat_obj$seurat_clusters %in% c(3, 5, 8, 9, 22, 27, 30) ~ "Myeloid",
    seurat_obj$seurat_clusters %in% c(7) ~ "Mast",
    seurat_obj$seurat_clusters %in% c(11) ~ "Fibroblast",
    seurat_obj$seurat_clusters %in% c(10) ~ "Endothelial",
    seurat_obj$seurat_clusters %in% c(15) ~ "Double_Pro_T_Myeloid",
    seurat_obj$seurat_clusters %in% c(28) ~ "Double_Mast_Epithelial",
    seurat_obj$seurat_clusters %in% c(29) ~ "Double_B_Epithelial",
    TRUE ~ "Unknown"
  )
  seurat_obj
}

#' Plot dot plot for marker genes.
#' @param seurat_obj Seurat object.
#' @param features Named list of marker genes.
#' @param group_by Column to group by.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_marker_dotplot <- function(seurat_obj, features, group_by, output_path,
                                width = 10, height = 4) {
  p <- DotPlot(object = seurat_obj, features = features, group.by = group_by) +
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
  pdf(file = output_path, width = width, height = height)
  print(p)
  dev.off()
  message("Saved marker dot plot to: ", output_path)
}

#' Plot UMAP of clusters and cell types.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_umap_annotation <- function(seurat_obj, output_path) {
  p3 <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE,
                group.by = "seurat_clusters") +
    theme(legend.position = "none")
  p4 <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE,
                group.by = "cell_type") +
    theme(legend.position = "none")
  pdf(file = output_path, width = 15, height = 6)
  print(p3 + p4)
  dev.off()
  message("Saved UMAP annotation to: ", output_path)
}

#' Plot final cell type UMAP with custom colors.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_final_umap <- function(seurat_obj, output_path) {
  custom_colors <- c(
    "#C6307C", "#D0AFC4", "#4991C1", "#89558D",
    "#AFC2D9", "#435B95", "#79B99D", "#F5A623"
  )
  p <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE,
               group.by = "cell_type") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  pdf(file = output_path, width = 7, height = 6)
  print(p)
  dev.off()
  message("Saved final UMAP to: ", output_path)
}

# ---- Main Function ----

#' Run cell annotation analysis.
#' @param input_path Path to filtered Seurat object.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return Final annotated Seurat object.
run_analysis <- function(
    input_path = file.path(RESULT_1_DIR, "results_data", "03_LUAD_merge_filter_decontx_double_gene.rds"),
    output_data_dir = file.path(RESULT_1_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_1_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  luad_merge <- load_filtered_data(input_path)
  luad_merge <- run_standard_workflow(luad_merge)

  # Unintegrated UMAP
  luad_merge <- RunUMAP(luad_merge, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

  # Harmony integration
  luad_merge <- run_harmony_integration(luad_merge, group_by = "Sample")
  plot_harmony_comparison(luad_merge, file.path(output_fig_dir, "01_harmony.pdf"))

  # Clustering
  luad_merge <- FindNeighbors(luad_merge, reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

  # Marker gene list
  marker_list <- list(
    "Epithelial" = c("EPCAM", "MUC1", "KRT19", "KRT18", "CDH1", "KRT8"),
    "T" = c("CD3D", "CD3E", "CD3G", "CD8A"),
    "NK" = c("NKG7", "KLRD1", "GNLY", "KLRB1"),
    "B" = c("CD79A", "MZB1", "JCHAIN", "MS4A1"),
    "Plasma" = c("IGHG1", "IGHA1"),
    "Myeloid" = c("CD68", "LYZ", "MARCO", "AIF1", "CTSD"),
    "Mast" = c("TPSB2", "MS4A2", "TPSAB1", "CPA3"),
    "Fibroblast" = c("COL1A1", "PDGFRA", "COL1A2", "ACTA2", "DCN", "LUM"),
    "Endothelial" = c("PECAM1", "VWF", "CDH5", "CLDN5", "RAMP2")
  )

  # Annotate cell types
  luad_merge <- annotate_cell_types(luad_merge)

  # Set cluster factor levels
  cluster_levels <- rev(c(
    2, 16, 17, 19, 20, 23, 0, 1, 6, 12, 13, 14, 18, 21, 26, 31,
    4, 24, 25, 3, 5, 8, 9, 22, 27, 30, 7, 11, 10, 15, 28, 29
  ))
  luad_merge$seurat_clusters <- factor(luad_merge$seurat_clusters, levels = cluster_levels)

  plot_marker_dotplot(luad_merge, marker_list, "seurat_clusters",
                      file.path(output_fig_dir, "02_major_cluster_marker.pdf"),
                      width = 12, height = 6)

  # Reorder cell types
  luad_merge$cell_type <- factor(luad_merge$cell_type,
    levels = rev(c(
      "Epithelial", "T/NK", "B", "Plasma", "Myeloid",
      "Mast", "Fibroblast", "Endothelial",
      "Double_Pro_T_Myeloid", "Double_Mast_Epithelial", "Double_B_Epithelial"
    ))
  )

  plot_marker_dotplot(luad_merge, marker_list, "cell_type",
                      file.path(output_fig_dir, "03_major_celltype_marker.pdf"),
                      width = 10, height = 4)

  plot_umap_annotation(luad_merge, file.path(output_fig_dir, "04_umap_clusters_cell_type.pdf"))

  saveRDS(luad_merge, file.path(output_data_dir, "05_LUAD_cell_anno.rds"))

  # Remove doublet clusters
  luad_merge <- subset(luad_merge, !cell_type %in% c(
    "Double_Pro_T_Myeloid", "Double_Mast_Epithelial", "Double_B_Epithelial"
  ))

  luad_merge$cell_type <- factor(luad_merge$cell_type,
    levels = rev(c("Epithelial", "T/NK", "B", "Plasma", "Myeloid",
                   "Mast", "Fibroblast", "Endothelial"))
  )

  plot_marker_dotplot(luad_merge, marker_list, "cell_type",
                      file.path(output_fig_dir, "05_celltype_marker.pdf"),
                      width = 10, height = 4)

  luad_merge$cell_type <- factor(luad_merge$cell_type,
    levels = c("Epithelial", "T/NK", "B", "Plasma", "Myeloid",
               "Mast", "Fibroblast", "Endothelial")
  )

  plot_final_umap(luad_merge, file.path(output_fig_dir, "06_umap_cell_type.pdf"))

  # Correct stage labels
  luad_merge$stage[which(luad_merge$stage == "Advanced")] <- "Advance"
  luad_merge$stage <- factor(luad_merge$stage, levels = c("Normal", "Early", "Advance"))

  saveRDS(luad_merge, file.path(output_data_dir, "06_LUAD_cell_anno_final.rds"))
  message("Cell annotation complete. Final object saved.")
  luad_merge
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in cell annotation: ", e$message)
      stop(e)
    }
  )
}
