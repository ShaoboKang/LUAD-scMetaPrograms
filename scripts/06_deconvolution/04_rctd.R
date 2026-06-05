# ============================================================================
# Script: 04_rctd
# Module: 06_deconvolution
# Description:
#   RCTD spatial deconvolution on Visium ST data using scRNA-seq reference.
#   Runs RCTD in full mode for 6 ST samples and saves results.
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

#' Run RCTD spatial deconvolution
#'
#' @param sc_path Path to Seurat RDS with MP assignments.
#' @param st_data_dir Directory with ST Seurat .rda files.
#' @param output_path Path to save RCTD results RDS.
#' @return List of RCTD objects.
run_rctd <- function(
    sc_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    st_data_dir = ST_DATA_DIR,
    output_path = file.path(RESULT_6_DIR, "results_data", "RCTD_visium_full.rds")
) {
  library(spacexr)
  library(Matrix)
  library(doParallel)
  library(ggplot2)

  message("Loading scRNA-seq reference...")
  luad_epi <- readRDS(sc_path)

  # Remap MP labels to consistent ordering
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4",
                 "MP_6", "MP_5", "MP_3", "MP_9")
  luad_epi$MP <- as.character(luad_epi$MP)
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))

  counts <- GetAssayData(luad_epi, assay = "RNA", slot = "counts")
  cell_types <- luad_epi$MP
  nUMI <- colSums(counts)

  reference <- Reference(counts, cell_types, nUMI)

  st_samples <- c("seurat_ST1", "seurat_ST2", "seurat_ST3",
                  "seurat_ST5", "seurat_ST6", "seurat_ST8")
  stages <- c("IAC", "IAC", "MIA", "AIS", "MIA", "AIS")

  RCTD <- list()
  for (i in seq_along(st_samples)) {
    message("Processing ", st_samples[i], " (", stages[i], ")...")
    load(file.path(st_data_dir, paste0(st_samples[i], ".rda")))

    st_counts <- seurat@assays$Spatial@counts
    coords <- seurat@images[["slice1"]]@coordinates[, c("row", "col")]
    colnames(coords) <- c("x", "y")

    st_nUMI <- colSums(st_counts)
    puck <- SpatialRNA(coords, st_counts, st_nUMI)

    barcodes <- colnames(puck@counts)
    plot_puck_continuous(
      puck, barcodes, puck@nUMI,
      ylimit = c(0, round(quantile(puck@nUMI, 0.9))),
      size = 2,
      title = paste0("nUMI: ", st_samples[i])
    )

    myRCTD <- create.RCTD(puck, reference, max_cores = 10)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
    RCTD[[paste0(st_samples[i], "_", stages[i])]] <- myRCTD
  }

  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(RCTD, output_path)
  message("[INFO] RCTD results saved to: ", output_path)
  invisible(RCTD)
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_rctd()
} else {
  message("[INFO] Script loaded. Call run_rctd() to execute.")
}
