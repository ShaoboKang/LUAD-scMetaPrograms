# ============================================================================
# Script: 01_merge_and_qc
# Module: 01_preprocessing
# Description:
#   Merge 6 scRNA-seq cohorts, run decontX ambient RNA removal,
#   scDblFinder doublet detection, mitochondrial/ribosome/HB filtering,
#   and generate a clean merged Seurat object.
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

#' Load all RDS files from a directory into a list.
#' @param dir_path Directory containing .rds files.
#' @return Named list of Seurat objects.
load_rds_files <- function(dir_path) {
  rds_files <- list.files(path = dir_path, pattern = "\\.rds$", full.names = FALSE)
  obj_list <- list()
  for (f in rds_files) {
    obj_list[[f]] <- readRDS(file.path(dir_path, f))
  }
  obj_list
}

#' Run decontX ambient RNA removal on a SingleCellExperiment.
#' @param sce SingleCellExperiment object.
#' @return Updated SCE with decontX results.
run_decontx <- function(sce) {
  sce <- decontX(sce, batch = sce$Sample)
  sce
}

#' Apply quality filters to a Seurat object.
#' @param seurat_obj Seurat object.
#' @return Filtered Seurat object.
apply_qc_filters <- function(seurat_obj) {
  seurat_obj <- PercentageFeatureSet(seurat_obj, "^MT-", col.name = "percent_mito")
  seurat_obj <- PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")
  seurat_obj <- PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent_mito < 10 & percent_hb < 1
  )
  seurat_obj
}

#' Detect and remove doublets using scDblFinder.
#' @param seurat_obj Seurat object.
#' @return Seurat object with doublets removed.
remove_doublets <- function(seurat_obj) {
  set.seed(123)
  sce_dbl <- as.SingleCellExperiment(seurat_obj)
  sce_dbl <- scDblFinder(sce_dbl, samples = "Sample", BPPARAM = MulticoreParam(10))
  seurat_obj$scDblFinder.class <- sce_dbl$scDblFinder.class
  seurat_obj <- subset(seurat_obj, scDblFinder.class == "singlet")
  seurat_obj
}

#' Remove mitochondrial, ribosomal, and hemoglobin genes.
#' @param seurat_obj Seurat object.
#' @return Seurat object with excluded genes removed.
remove_excluded_genes <- function(seurat_obj) {
  gene_set <- rownames(seurat_obj)
  mito_genes <- gene_set[grep("^MT-", gene_set, ignore.case = TRUE)]
  ribo_genes <- gene_set[grep("^RP[SL]", gene_set, ignore.case = TRUE)]
  hb_genes <- gene_set[grep("^HB[^(P)]", gene_set, ignore.case = TRUE)]
  genes_to_keep <- setdiff(gene_set, c(mito_genes, ribo_genes, hb_genes))
  seurat_obj <- seurat_obj[genes_to_keep, ]
  seurat_obj
}

#' Remove genes expressed in fewer than a minimum number of cells.
#' @param seurat_obj Seurat object.
#' @param min_cells Minimum number of cells with expression > 0.
#' @return Filtered Seurat object.
filter_low_expression_genes <- function(seurat_obj, min_cells = 10) {
  counts_layer <- seurat_obj@assays[["RNA"]]@layers[["counts"]]
  expressed_cells <- rowSums(counts_layer > 0)
  features_to_keep <- rownames(seurat_obj)[expressed_cells > min_cells]
  subset(seurat_obj, features = features_to_keep)
}

# ---- Main Function ----

#' Run the merge and QC pipeline.
#' @param data_dir Directory containing input .rds files.
#' @param output_dir Directory for output files.
#' @return The final filtered Seurat object.
run_analysis <- function(data_dir = SC_DATA_DIR,
                         output_dir = file.path(RESULT_1_DIR, "results_data")) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Load and merge all cohorts
  message("Loading RDS files from ", data_dir)
  luad_merge_list <- load_rds_files(data_dir)
  message("Merging ", length(luad_merge_list), " cohorts...")
  luad_merge <- merge(luad_merge_list[[1]], luad_merge_list[2:length(luad_merge_list)]) %>%
    JoinLayers()
  message("Merged object: ", ncol(luad_merge), " cells")

  # Convert to SingleCellExperiment and run decontX
  message("Running decontX for ambient RNA removal...")
  sce <- as.SingleCellExperiment(luad_merge)
  sce <- run_decontx(sce)

  # Update counts and contamination info
  luad_merge <- SetAssayData(
    object = luad_merge,
    assay = "RNA",
    layer = "counts",
    new.data = round(decontXcounts(sce))
  )
  luad_merge$decontX_contamination <- sce$decontX_contamination

  # Filter high contamination cells
  luad_merge <- subset(luad_merge, subset = decontX_contamination < 0.2)
  message("After decontX filtering: ", ncol(luad_merge), " cells")

  # Apply QC filters
  message("Applying QC filters...")
  luad_merge <- apply_qc_filters(luad_merge)
  message("After QC filtering: ", ncol(luad_merge), " cells")

  # Remove doublets
  message("Running scDblFinder...")
  luad_merge <- remove_doublets(luad_merge)
  message("After doublet removal: ", ncol(luad_merge), " cells")

  # Save doublet info
  saveRDS(luad_merge$scDblFinder.class,
          file = file.path(output_dir, "02_scDblFinder.rds"))

  # Remove excluded gene sets
  message("Removing mitochondrial, ribosomal, and hemoglobin genes...")
  luad_merge <- remove_excluded_genes(luad_merge)

  # Remove low-expression genes
  message("Removing genes expressed in <= 10 cells...")
  luad_merge <- filter_low_expression_genes(luad_merge, min_cells = 10)

  # Save final object
  out_file <- file.path(output_dir, "03_LUAD_merge_filter_decontx_double_gene.rds")
  saveRDS(luad_merge, out_file)
  message("Final object saved to: ", out_file)

  luad_merge
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in merge and QC pipeline: ", e$message)
      stop(e)
    }
  )
}
