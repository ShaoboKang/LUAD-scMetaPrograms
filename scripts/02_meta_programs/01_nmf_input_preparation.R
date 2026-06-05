# ============================================================================
# Script: 01_nmf_input_preparation
# Module: 02_meta_programs
# Description:
#   Extract tumor epithelial cells, filter samples with >100 cells,
#   run standard Seurat workflow, and save NMF input object.
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

# ---- Helper Functions ----

#' Load final epithelial annotation.
#' @param path Path to .rds file.
#' @return Seurat object.
load_epithelial_data <- function(path) {
  message("Loading epithelial data from: ", path)
  readRDS(path)
}

#' Subset to tumor epithelial cells and filter small samples.
#' @param seurat_obj Seurat object.
#' @param min_cells Minimum cells per sample. Default 100.
#' @return Filtered Seurat object.
filter_tumor_samples <- function(seurat_obj, min_cells = 100) {
  seurat_obj <- subset(seurat_obj, Tissue == "Tumor")
  seurat_obj <- subset(seurat_obj, Epi_subtype == "Tumor")
  sample_cell_count <- table(seurat_obj$Sample)
  large_samples <- names(sample_cell_count[sample_cell_count > min_cells])
  seurat_obj <- subset(seurat_obj, subset = Sample %in% large_samples)
  message("Retained ", length(large_samples), " samples with >", min_cells, " cells")
  seurat_obj
}

#' Run standard Seurat preprocessing.
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

# ---- Main Function ----

#' Prepare NMF input data.
#' @param input_path Path to final epithelial annotation.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return Prepared Seurat object.
run_analysis <- function(
    input_path = file.path(RESULT_1_DIR, "results_data", "08_LUAD_epi_anno_final.rds"),
    output_data_dir = file.path(RESULT_2_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  luad_epi <- load_epithelial_data(input_path)
  luad_epi <- filter_tumor_samples(luad_epi, min_cells = 100)
  luad_epi <- run_standard_workflow(luad_epi)

  luad_epi <- RunUMAP(luad_epi, dims = 1:30, reduction = "pca")

  pdf(file = file.path(output_fig_dir, "01_umap_tumor_sample.pdf"), width = 7, height = 6)
  print(DimPlot(luad_epi, reduction = "umap", label = TRUE, group.by = "Sample", repel = TRUE) +
    theme(legend.position = "none"))
  dev.off()
  message("Saved tumor sample UMAP")

  saveRDS(luad_epi, file = file.path(output_data_dir, "01_LUAD_epi_anno_final_NMF_input.rds"))
  message("NMF input preparation complete.")
  luad_epi
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in NMF input preparation: ", e$message)
      stop(e)
    }
  )
}
