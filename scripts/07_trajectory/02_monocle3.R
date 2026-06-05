# ============================================================================
# Script: 02_monocle3
# Module: 07_trajectory
# Description: Monocle3 trajectory inference and TF activity along pseudotime.
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

#' Import integrated UMAP coordinates from Seurat into Monocle3 CDS
#'
#' @param cds Monocle3 cell_data_set object.
#' @param seurat_obj Seurat object with umap.harmony reduction.
#' @return CDS with updated UMAP coordinates.
import_seurat_umap <- function(cds, seurat_obj) {
  cds@int_colData$reducedDims$UMAP <- Embeddings(seurat_obj, reduction = "umap.harmony")
  return(cds)
}

# ---- Main Function ----

#' Run Monocle3 trajectory and TF activity analysis
#'
#' Builds a Monocle3 object from Seurat counts, imports UMAP coordinates,
#' learns the trajectory graph, orders cells by pseudotime, and links
#' SCENIC TF activities to pseudotime.
#'
#' @param sc_path Path to Seurat RDS.
#' @param loom_path Path to SCENIC loom file.
#' @param output_dir Directory for outputs.
#' @return Ordered Monocle3 CDS object.
run_analysis <- function(
    sc_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds"),
    loom_path = file.path(RESULT_4_DIR, "scenic", "output", "out_SCENIC.loom"),
    output_dir = file.path(RESULT_7_DIR, "results_data")
) {
  library(igraph)
  library(monocle3)
  library(Seurat)
  library(ggplot2)
  library(SCopeLoomR)
  library(lemon)
  library(patchwork)

  message("Loading Seurat object...")
  luad_epi <- readRDS(sc_path)

  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))

  # Construct monocle3 CDS
  cds <- new_cell_data_set(
    expression_data = GetAssayData(luad_epi, slot = "counts"),
    cell_metadata = luad_epi@meta.data,
    gene_metadata = data.frame(
      gene_short_name = row.names(luad_epi),
      row.names = row.names(luad_epi)
    )
  )

  # Preprocessing: normalization, batch correction, and dimensionality reduction
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- align_cds(cds, alignment_group = "Sample")
  cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP")

  # Import integrated UMAP from Seurat
  cds <- import_seurat_umap(cds, luad_epi)
  plot_cells(cds, reduction_method = "UMAP", color_cells_by = "MP")

  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)

  root_cell <- colnames(cds)[which(cds@clusters@listData[["UMAP"]][["clusters"]] == c("2"))]
  cds <- order_cells(cds, root_cells = root_cell)

  fig_dir <- file.path(RESULT_7_DIR, "results_figure")
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  pdf(file.path(fig_dir, "01_monocle3.pdf"), width = 5, height = 4)
  print(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   show_trajectory_graph = TRUE, label_roots = FALSE,
                   label_leaves = FALSE, label_branch_points = TRUE) +
          scale_color_viridis_c(option = "D"))
  dev.off()

  pdf(file.path(fig_dir, "02_monocle3_stage.pdf"), width = 6, height = 4)
  print(plot_cells(cds, color_cells_by = "stage", label_cell_groups = FALSE,
                   show_trajectory_graph = TRUE, label_roots = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
          scale_color_manual(values = c("#1F77B4", "#D62728")))
  dev.off()

  pdf(file.path(fig_dir, "03_monocle3_cluster.pdf"), width = 6, height = 4)
  print(plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE,
                   show_trajectory_graph = TRUE, label_roots = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
          scale_color_manual(values = c(
            "#E64B35FF", "#4DBBD5FF", "#91D1C2FF", "#00A087FF",
            "#3C5488FF", "#F39B7FFF", "grey"
          )))
  dev.off()

  plot_cells(cds, color_cells_by = "MP", label_cell_groups = FALSE,
             show_trajectory_graph = TRUE, label_roots = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE) +
    scale_color_manual(values = c(
      "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
      "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
    ))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cds, file = file.path(output_dir, "03_cds_monocle3.rds"))

  # TF activity along trajectory
  message("Loading SCENIC loom and extracting TF activities...")
  loom <- open_loom(loom_path)
  regulons_incid_mat <- get_regulons(loom, column.attr.name = "Regulons")
  regulon_auc <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  tf_activity <- regulon_auc@assays@data@listData[["AUC"]]

  tf_activate_pseudotime <- data.frame(
    cell_id = colnames(cds),
    pseudotime = pseudotime(cds),
    branch = cds@clusters@listData[["UMAP"]][["clusters"]]
  )
  tf_activate_pseudotime <- cbind(tf_activate_pseudotime, t(tf_activity))

  saveRDS(tf_activate_pseudotime, file = file.path(output_dir, "TF_activate_pseudotime.rds"))

  tf_auc <- t(tf_activity[c("XBP1(+)", "JUN(+)", "PPARG(+)", "FOSL2(+)", "FOXM1(+)", "VPS72(+)"), ])
  luad_epi <- AddMetaData(luad_epi, tf_auc,
                          col.name = paste0(gsub("\\(.*\\)", "", colnames(tf_auc)), " Activity"))

  features <- c("XBP1 Activity", "JUN Activity", "PPARG Activity",
               "FOSL2 Activity", "FOXM1 Activity", "VPS72 Activity")

  plots <- lapply(features, function(feature) {
    FeaturePlot(luad_epi, reduction = "umap.harmony", features = feature) +
      scale_color_viridis_c(option = "D") +
      ggtitle(feature)
  })

  pdf(file.path(fig_dir, "04_TF_auc_pesudo.pdf"), width = 15, height = 8)
  print(wrap_plots(plots, nrow = 2, ncol = 3))
  dev.off()

  message("Monocle3 analysis complete.")
  return(cds)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
