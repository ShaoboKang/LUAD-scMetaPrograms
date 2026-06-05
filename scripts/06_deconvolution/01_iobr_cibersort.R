# ============================================================================
# Script: 01_iobr_cibersort
# Module: 06_deconvolution
# Description: Run CIBERSORT deconvolution on bulk RNA-seq using MP signatures.
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

#' Reorder MP labels in a Seurat object
#'
#' @param seurat_obj A Seurat object containing MP labels.
#' @param old_order Character vector defining the desired MP order.
#' @return The Seurat object with reordered MP labels.
reorder_mp_labels <- function(seurat_obj, old_order) {
  seurat_obj$MP <- as.character(seurat_obj$MP)
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  seurat_obj$MP <- unname(mp_map[seurat_obj$MP])
  seurat_obj$MP <- factor(seurat_obj$MP, levels = paste0("MP_", seq_along(old_order)))
  return(seurat_obj)
}

#' Generate a CIBERSORT reference matrix from a Seurat object
#'
#' @param seurat_obj A Seurat object with normalized data and MP labels.
#' @return A reference matrix for CIBERSORT.
generate_cibersort_reference <- function(seurat_obj) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  reference <- IOBR::generateRef_seurat(
    sce = seurat_obj,
    celltype = "MP",
    slot_out = "data"
  )
  return(reference)
}

#' Run CIBERSORT deconvolution and attach clinical group labels
#'
#' @param mrna_exp Bulk mRNA expression matrix.
#' @param reference CIBERSORT reference matrix.
#' @param clinical Clinical data frame containing stage information.
#' @return Data frame of CIBERSORT results.
run_cibersort <- function(mrna_exp, reference, clinical) {
  result <- IOBR::deconvo_tme(
    eset = mrna_exp,
    reference = reference,
    method = "cibersort",
    arrays = FALSE,
    absolute.mode = FALSE,
    perm = 100
  )
  result <- as.data.frame(result)
  rownames(result) <- result$ID
  result$group <- clinical$stage
  result$group <- factor(result$group, levels = c("normal", "early", "advance"))
  colnames(result)[2:10] <- paste0("MP_", 1:9)
  result[, 2:10] <- apply(result[, 2:10], 2, as.numeric)
  return(result)
}

#' Create boxplots comparing MP fractions across clinical stages
#'
#' @param cibersort_df Data frame of CIBERSORT results.
#' @return A patchwork object of boxplots.
plot_cibersort_boxplots <- function(cibersort_df) {
  target_columns <- setdiff(
    colnames(cibersort_df),
    c("ID", "group", "P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT")
  )

  boxplot_list <- list()
  for (i in target_columns) {
    boxplot_list[[i]] <- ggpubr::ggboxplot(
      cibersort_df,
      x = "group",
      y = i,
      fill = "group",
      outlier.shape = NA,
      title = i
    ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank()
      ) +
      ggpubr::stat_compare_means(
        comparisons = list(c("early", "normal"), c("advance", "early"), c("advance", "normal")),
        method = "wilcox.test"
      ) +
      ggplot2::scale_fill_manual(
        values = c("normal" = "#fce38a", "early" = "#1F77B4", "advance" = "#D62728")
      )
  }
  patchwork::wrap_plots(boxplot_list, ncol = 3)
}

# ---- Main Function ----

#' Run CIBERSORT deconvolution analysis and save outputs
#'
#' @param mp_list_path Path to MP list RData file.
#' @param seurat_path Path to assigned MP Seurat RDS file.
#' @param bulk_path Path to bulk clinical RData file.
#' @param output_dir Output directory for results.
run_analysis <- function(mp_list_path, seurat_path, bulk_path, output_dir) {
  message("Loading MP list from: ", mp_list_path)
  load(mp_list_path)

  message("Loading clinical data from: ", bulk_path)
  load(bulk_path)

  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  MP_list <- MP_list[old_order]
  names(MP_list) <- paste0("MP_", seq_along(MP_list))

  message("Loading single-cell data from: ", seurat_path)
  LUAD_Epi_assigned_MP <- readRDS(seurat_path)
  LUAD_Epi_assigned_MP <- reorder_mp_labels(LUAD_Epi_assigned_MP, old_order)

  message("Generating CIBERSORT reference")
  lm22 <- generate_cibersort_reference(LUAD_Epi_assigned_MP)

  message("Running CIBERSORT deconvolution")
  cibersort <- run_cibersort(mRNA_exp, lm22, clinical)

  message("Creating boxplots")
  boxplot_figure <- plot_cibersort_boxplots(cibersort)

  figure_path <- file.path(output_dir, "results_figure", "01_cibersort.pdf")
  message("Saving figure to: ", figure_path)
  pdf(file = figure_path, width = 12, height = 12)
  print(boxplot_figure)
  dev.off()

  data_path <- file.path(output_dir, "results_data", "01_cibersort.rds")
  message("Saving results to: ", data_path)
  saveRDS(cibersort, file = data_path)

  message("CIBERSORT analysis complete")
}

# ---- Entry Point ----
if (!interactive()) {
  run_analysis(
    mp_list_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    seurat_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA", "TCGA_sur_data_stage.Rdata"),
    output_dir = RESULT_6_DIR
  )
}
