# ============================================================================
# Script: 06_gene_switch_analysis
# Utility: geneswitch.R
# Description: GeneSwitch analysis along Monocle3 lineage branches to
#   identify gene switching events.
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

#' Convert Seurat object to GeneSwitches SingleCellExperiment
#'
#' @param seurat_rds Seurat object.
#' @param clusters Vector of cluster IDs to subset.
#' @param tf_names TF gene vector.
#' @param mp_gene MP gene vector.
#' @return SingleCellExperiment object.
seurat_to_gene_switch_sce <- function(seurat_rds, clusters, tf_names, mp_gene) {
  seurat_rds <- subset(seurat_rds, mono_cluster %in% clusters)
  sce <- SingleCellExperiment(
    assays = List(expdata = GetAssayData(seurat_rds, layer = "data"))
  )
  sce <- sce[unique(c(tf_names, mp_gene)), ]
  colData(sce)$Pseudotime <- seurat_rds@meta.data$pseudotime
  umap <- seurat_rds@reductions$umap.unintegrated@cell.embeddings
  reducedDims(sce) <- SimpleList(UMAP = umap)
  return(sce)
}

#' Run GeneSwitch analysis for both lineages
#'
#' @param mp_list_path Path to MP top50 RData.
#' @param tf_path Path to top5 TF RDS.
#' @param cds_path Path to Monocle3 CDS RDS.
#' @param sc_path Path to Seurat RDS.
#' @param binarize_path Path to binarize_exp utility.
#' @param output_fig_dir Directory for figures.
#' @return List of GeneSwitch plots.
run_geneswitch_analysis <- function(
    mp_list_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    tf_path = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"),
    cds_path = file.path(RESULT_7_DIR, "results_data", "03_cds_monocle3.rds"),
    sc_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds"),
    binarize_path = file.path(CODE_DIR, "scripts", "07_trajectory", "utils", "binarize_exp.R"),
    output_fig_dir = file.path(RESULT_7_DIR, "results_figure")
) {
  if (exists("memory.limit", mode = "function")) {
    memory.limit(size = 900000000)
  }
  library(GeneSwitches)
  library(SingleCellExperiment)
  library(ggplot2)
  library(Seurat)
  library(patchwork)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading data...")
  load(mp_list_path)
  tf_names <- readRDS(tf_path)
  tf_names <- unique(gsub("\\(.*\\)", "", tf_names$Topic))
  mp_gene <- unique(unlist(MP_list))

  cds <- readRDS(cds_path)
  luad_epi <- readRDS(sc_path)
  luad_epi@meta.data$pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
  luad_epi@meta.data$mono_cluster <- cds@clusters@listData[["UMAP"]][["clusters"]]

  source(binarize_path)

  # Lineage 1: clusters 2,4,1
  message("Processing lineage 1 (clusters 2,4,1)...")
  sce_p1 <- seurat_to_gene_switch_sce(luad_epi, c(2, 4, 1), tf_names, mp_gene)
  sce_p1 <- binarize_exp(sce_p1, ncores = 1, fix_cutoff = TRUE, binarize_cutoff = 0.2)
  sce_p1 <- find_switch_logistic_fastglm(sce_p1, show_warning = FALSE)

  sg_tf <- filter_switchgenes(sce_p1, allgenes = TRUE)
  sg_tf <- sg_tf[intersect(
    sg_tf@rownames,
    c(unique(unlist(MP_list[c("MP_2", "MP_3", "MP_7", "MP_9")])), tf_names)
  ), ]
  sg_tf$feature_type <- "MP_gene"
  sg_tf$feature_type[match(intersect(sg_tf@rownames, tf_names), sg_tf@rownames)] <- "TFs"

  p1 <- plot_timeline_ggplot(sg_tf, timedata = sce_p1$Pseudotime,
                             txtsize = 3, color_by = "feature_type") +
    scale_color_manual(values = c("#1F77B4", "#D62728"))

  # Lineage 2: clusters 2,6,3
  message("Processing lineage 2 (clusters 2,6,3)...")
  sce_p2 <- seurat_to_gene_switch_sce(luad_epi, c(2, 6, 3), tf_names, mp_gene)
  sce_p2 <- binarize_exp(sce_p2, ncores = 1, fix_cutoff = TRUE, binarize_cutoff = 0.2)
  sce_p2 <- find_switch_logistic_fastglm(sce_p2, show_warning = FALSE)

  sg_tf2 <- filter_switchgenes(sce_p2, allgenes = TRUE)
  sg_tf2 <- sg_tf2[intersect(
    sg_tf2@rownames,
    c(unique(unlist(MP_list[c("MP_2", "MP_7", "MP_5", "MP_1", "MP_8")])), tf_names)
  ), ]
  sg_tf2$feature_type <- "MP_gene"
  sg_tf2$feature_type[match(intersect(sg_tf2@rownames, tf_names), sg_tf2@rownames)] <- "TFs"

  p2 <- plot_timeline_ggplot(sg_tf2, timedata = sce_p2$Pseudotime,
                             txtsize = 3, color_by = "feature_type") +
    scale_color_manual(values = c("#1F77B4", "#D62728"))

  pdf(file.path(output_fig_dir, "geneswitch.pdf"), width = 7, height = 8)
  print(p1 / p2)
  dev.off()

  message("GeneSwitch analysis complete.")
  return(list(lineage1 = p1, lineage2 = p2))
}
# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_geneswitch_analysis()
}
