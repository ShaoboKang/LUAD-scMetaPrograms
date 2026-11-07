# Restore the preparation step recorded in Result_8/code/02_scissor.R.
# Keep legacy MP labels here; downstream analyses apply the state mapping once.
source("config.R")
suppressPackageStartupMessages({library(Seurat); library(harmony)})

#' Prepare the shared integrated MP object for trajectory and Scissor analyses.
#' @param input_path Assigned malignant-cell object with legacy MP labels.
#' @param output_path Output RDS containing PCA and Harmony UMAP embeddings.
#' @return Integrated Seurat object, invisibly.
run_analysis <- function(
    input_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    output_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds")) {
  obj <- readRDS(input_path)
  stopifnot("Sample" %in% names(obj@meta.data), "MP" %in% names(obj@meta.data))
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  obj <- RunHarmony(obj, group.by.vars = "Sample")
  obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, output_path)
  invisible(obj)
}
if (!interactive()) run_analysis()
