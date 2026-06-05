# ============================================================================
# Script: 04_epithelial_subtyping
# Module: 01_preprocessing
# Description:
#   Subset epithelial cells, re-cluster with Harmony, define epithelial
#   subtypes (AT2, AT1, Club, Ciliated, Tumor), and remove doublet clusters.
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
  library(ggsci)

# ---- Helper Functions ----

#' Load annotated data and subset to epithelial cells.
#' @param path Path to annotated Seurat object.
#' @return Epithelial subset Seurat object.
load_epithelial_cells <- function(path) {
  message("Loading data from: ", path)
  readRDS(path) %>% subset(cell_type == "Epithelial")
}

#' Run standard Seurat workflow.
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

#' Integrate epithelial cells with Harmony.
#' @param seurat_obj Seurat object with PCA.
#' @param group_by Column to integrate by. Default "Study".
#' @return Harmony-integrated Seurat object.
run_harmony_integration <- function(seurat_obj, group_by = "Study") {
  set.seed(123)
  seurat_obj <- RunHarmony(seurat_obj, group.by = group_by)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
  seurat_obj
}

#' Plot UMAP before and after integration.
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

#' Annotate epithelial subtypes from cluster assignments.
#' @param seurat_obj Seurat object with seurat_clusters.
#' @return Seurat object with Epi_subtype column.
annotate_epithelial_subtypes <- function(seurat_obj) {
  seurat_obj$Epi_subtype <- case_when(
    seurat_obj$seurat_clusters %in% c(1) ~ "AT2",
    seurat_obj$seurat_clusters %in% c(10) ~ "AT1",
    seurat_obj$seurat_clusters %in% c(6, 22) ~ "Club",
    seurat_obj$seurat_clusters %in% c(7) ~ "Ciliated",
    seurat_obj$seurat_clusters %in% c(11, 26) ~ "Double_T/NK",
    TRUE ~ "Tumor"
  )
  seurat_obj
}

#' Plot dot plot for epithelial markers.
#' @param seurat_obj Seurat object.
#' @param features Named list of marker genes.
#' @param group_by Column to group by.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_marker_dotplot <- function(seurat_obj, features, group_by, output_path,
                                width = 10, height = 6) {
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

#' Plot UMAP of clusters and epithelial subtypes.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_umap_subtypes <- function(seurat_obj, output_path) {
  custom_colors <- rev(c(
    "#F39B7FFF", "#4DBBD5FF", "#00A087FF",
    "#3C5488FF", "#E64B35FF", "grey"
  ))
  p3 <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE,
                group.by = "seurat_clusters") +
    theme(legend.position = "none")
  p4 <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE,
                group.by = "Epi_subtype") +
    scale_color_manual(values = custom_colors) +
    theme(legend.position = "none")
  pdf(file = output_path, width = 15, height = 6)
  print(p3 + p4)
  dev.off()
  message("Saved UMAP subtypes to: ", output_path)
}

#' Plot final epithelial subtype UMAP.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_final_umap <- function(seurat_obj, output_path) {
  custom_colors <- rev(c(
    "#F39B7FFF", "#4DBBD5FF", "#00A087FF",
    "#3C5488FF", "#E64B35FF"
  ))
  p <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE,
               group.by = "Epi_subtype") +
    scale_color_manual(values = custom_colors) +
    theme(legend.position = "none")
  pdf(file = output_path, width = 8, height = 6)
  print(p)
  dev.off()
  message("Saved final UMAP to: ", output_path)
}

#' Plot barplot of tissue composition per cluster.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
plot_tissue_barplot <- function(seurat_obj, output_path) {
  cellratio <- as.data.frame(prop.table(table(seurat_obj$Tissue, seurat_obj$seurat_clusters), margin = 2))
  cellratio$Var1 <- factor(cellratio$Var1, levels = c("Tumor", "Normal"))
  p <- ggplot(cellratio) +
    geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
    theme_classic() +
    scale_color_jco() +
    scale_fill_jco() +
    labs(x = "Cluster", y = "Ratio") +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
  pdf(file = output_path, width = 10, height = 4)
  print(p)
  dev.off()
  message("Saved tissue barplot to: ", output_path)
}

# ---- Main Function ----

#' Run epithelial subtyping analysis.
#' @param input_path Path to annotated Seurat object.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return Final epithelial Seurat object.
run_analysis <- function(
    input_path = file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds"),
    output_data_dir = file.path(RESULT_1_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_1_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  luad_epi <- load_epithelial_cells(input_path)
  luad_epi <- run_standard_workflow(luad_epi)

  # Unintegrated UMAP
  luad_epi <- RunUMAP(luad_epi, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

  # Harmony integration
  luad_epi <- run_harmony_integration(luad_epi, group_by = "Study")
  plot_harmony_comparison(luad_epi, file.path(output_fig_dir, "09_epi_harmony.pdf"))

  # Clustering
  luad_epi <- FindNeighbors(luad_epi, reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

  # Annotate subtypes
  luad_epi <- annotate_epithelial_subtypes(luad_epi)

  marker_list <- list(
    "AT2" = c("SFTPB", "SFTPC", "SFTPD", "SFTPA1"),
    "AT1" = c("AGER", "PDPN", "CLIC5", "CAV1"),
    "Club" = c("SCGB3A2", "SCGB3A1", "SCGB1A1", "MUC5B"),
    "Ciliated" = c("TUBA1A", "FOXJ1", "CAPS", "CCDC78"),
    "Tumor" = c("MDK", "TIMP1"),
    "T/NK" = c("CD3D", "CD3E", "CD3G", "NKG7", "GNLY")
  )

  # Set cluster levels
  cluster_levels <- rev(c(
    1, 10, 6, 22, 7, 0, 2, 3, 4, 5, 8, 9,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 27, 11, 26
  ))
  luad_epi$seurat_clusters <- factor(luad_epi$seurat_clusters, levels = cluster_levels)

  plot_marker_dotplot(luad_epi, marker_list, "seurat_clusters",
                      file.path(output_fig_dir, "10_epi_cluster_marker.pdf"),
                      width = 10, height = 6)

  luad_epi$Epi_subtype <- factor(luad_epi$Epi_subtype,
    levels = rev(c("AT2", "AT1", "Club", "Ciliated", "Tumor", "Double_T/NK"))
  )

  plot_marker_dotplot(luad_epi, marker_list, "Epi_subtype",
                      file.path(output_fig_dir, "11_epi_subtype_marker.pdf"),
                      width = 8, height = 4)

  plot_umap_subtypes(luad_epi, file.path(output_fig_dir, "12_umap_clusters_cell_type.pdf"))

  saveRDS(luad_epi, file.path(output_data_dir, "07_LUAD_epi_anno.rds"))

  # Remove doublet subtype
  luad_epi <- subset(luad_epi, !Epi_subtype %in% "Double_T/NK")

  plot_final_umap(luad_epi, file.path(output_fig_dir, "13_umap_celltype_final.pdf"))

  luad_epi$seurat_clusters <- luad_epi$RNA_snn_res.0.5
  luad_epi$seurat_clusters <- factor(luad_epi$seurat_clusters, levels = c(0:10, 12:25, 27))

  saveRDS(luad_epi, file.path(output_data_dir, "08_LUAD_epi_anno_final.rds"))

  plot_tissue_barplot(luad_epi, file.path(output_fig_dir, "14_barplot_cluster_ratio_final.pdf"))

  message("Epithelial subtyping complete.")
  luad_epi
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in epithelial subtyping: ", e$message)
      stop(e)
    }
  )
}
