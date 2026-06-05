# ============================================================================
# Script: 08_mp9_emt_stemness_correlation
# Module: 02_meta_programs
# Description:
#   Compute overlap between MP gene sets and published MP41 hallmark gene sets
#   and visualize as a heatmap.
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
  library(xlsx)

# ---- Helper Functions ----

#' Load MP gene lists.
#' @param rdata_path Path to RData file.
#' @return Data frame of MP genes.
load_mp_list <- function(rdata_path) {
  message("Loading MP list from: ", rdata_path)
  load(rdata_path)
  new_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_list <- MP_list[new_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))
  as.data.frame(mp_list)
}

#' Read MP41 hallmark gene sets from Excel.
#' @param xlsx_path Path to Excel file.
#' @return Data frame of gene sets.
load_mp41_hallmarks <- function(xlsx_path) {
  message("Loading MP41 hallmarks from: ", xlsx_path)
  read.xlsx(xlsx_path, sheetIndex = 1)
}

#' Compute pairwise intersection matrix.
#' @param mp_list Data frame of MP genes.
#' @param hallmark_df Data frame of hallmark gene sets.
#' @return Data frame of intersection sizes.
compute_intersection_matrix <- function(mp_list, hallmark_df) {
  result_matrix <- matrix(0, nrow = ncol(mp_list), ncol = ncol(hallmark_df))
  for (i in 1:ncol(mp_list)) {
    for (j in 1:ncol(hallmark_df)) {
      mp_genes <- mp_list[[i]]
      data_genes <- hallmark_df[[j]]
      result_matrix[i, j] <- length(intersect(mp_genes, data_genes))
    }
  }
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- colnames(hallmark_df)
  rownames(result_df) <- colnames(mp_list)
  result_df
}

#' Plot intersection heatmap.
#' @param result_df Data frame of intersections.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_intersection_heatmap <- function(result_df, output_path, width = 14, height = 6) {
  bk <- c(seq(-4, -0.1, by = 0.01), seq(0, 4, by = 0.01))
  pdf(file = output_path, width = width, height = height)
  pheatmap::pheatmap(
    result_df,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 10,
    use_raster = FALSE,
    border_color = "black",
    color = c(
      colorRampPalette(colors = c("white", "#091A8C"))(length(bk) / 2),
      colorRampPalette(colors = c("#091A8C", "#B30000"))(length(bk) / 2)
    )
  )
  dev.off()
  message("Saved intersection heatmap to: ", output_path)
}

# ---- Main Function ----

#' Run MP9 EMT/stemness correlation analysis.
#' @param rdata_path Path to MP RData file.
#' @param xlsx_path Path to MP41 hallmark Excel file.
#' @param output_fig_dir Directory for output figures.
#' @return Intersection data frame.
run_analysis <- function(
    rdata_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    xlsx_path = file.path(MSIGDB_DIR, "MP41_hallmarker.xlsx"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_list <- load_mp_list(rdata_path)
  hallmark_df <- load_mp41_hallmarks(xlsx_path)

  result_df <- compute_intersection_matrix(mp_list, hallmark_df)

  plot_intersection_heatmap(result_df, file.path(output_fig_dir, "17_MP1-9_cor_MP41.pdf"))

  message("MP9 EMT/stemness correlation analysis complete.")
  invisible(result_df)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP9 correlation analysis: ", e$message)
      stop(e)
    }
  )
}
