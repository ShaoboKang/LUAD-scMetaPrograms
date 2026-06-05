# ============================================================================
# Script: 03_tumor_ratio_clustering
# Module: 06_deconvolution
# Description: Cluster tumor samples by MP proportions and MP scores.
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

#' Compute MP proportions per sample
#'
#' @param seurat_obj A Seurat object with MP labels and Sample metadata.
#' @return Data frame with MP, sample, count, and proportion columns.
compute_mp_proportions <- function(seurat_obj) {
  count_table <- as.data.frame(table(seurat_obj$MP, seurat_obj$Sample))
  names(count_table) <- c("MP", "sample", "count")
  count_table <- count_table %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
  return(count_table)
}

#' Build sample annotation data frame from Seurat metadata
#'
#' @param seurat_obj A Seurat object with stage and AJCC_stage metadata.
#' @param sample_names Character vector of sample names matching row names.
#' @return Data frame with stage and ajcc_stage annotations.
build_sample_annotation <- function(seurat_obj, sample_names) {
  annotation_col <- data.frame(
    row.names = sample_names,
    stage = seurat_obj$stage[match(sample_names, seurat_obj$Sample)],
    ajcc_stage = seurat_obj$AJCC_stage[match(sample_names, seurat_obj$Sample)]
  )
  annotation_col$ajcc_stage <- dplyr::case_when(
    annotation_col$ajcc_stage %in% c("MIA", "AIS", "IA", "IA3", "IAC", "IB", "IB-IIA") ~ "I",
    annotation_col$ajcc_stage %in% c("IIA", "IIB") ~ "II",
    annotation_col$ajcc_stage %in% c("IIIA", "IIIB", "IIIB ") ~ "III",
    annotation_col$ajcc_stage %in% c("IV") ~ "IV"
  )
  return(annotation_col)
}

#' Plot sample-by-MP proportion heatmap
#'
#' @param proportion_wide Wide-format proportion matrix.
#' @param annotation_col Sample annotation data frame.
#' @param ann_colors List of annotation color mappings.
#' @param output_path Path to save PDF figure.
plot_proportion_heatmap <- function(proportion_wide, annotation_col, ann_colors, output_path) {
  pdf(output_path, width = 12, height = 4)
  pheatmap::pheatmap(
    t(proportion_wide),
    clustering_method = "ward.D2",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = "Tumor cell state ratio across sc-samples"
  )
  dev.off()
}

#' Plot sample-by-MP score heatmap sorted by AJCC stage
#'
#' @param mp_score Data frame of per-sample MP scores.
#' @param annotation_col Sample annotation data frame.
#' @param ann_colors List of annotation color mappings.
#' @param output_path Path to save PDF figure.
plot_score_heatmap <- function(mp_score, annotation_col, ann_colors, output_path) {
  annotation_col <- annotation_col[order(annotation_col$ajcc_stage), , drop = FALSE]
  mp_score <- mp_score[rownames(annotation_col), , drop = FALSE]
  mp_score$ajcc_stage <- annotation_col$ajcc_stage

  mp_score_sorted <- mp_score %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::mutate(ajcc_stage = factor(ajcc_stage, levels = c("I", "II", "III", "IV"))) %>%
    dplyr::group_by(ajcc_stage) %>%
    dplyr::arrange(MP_1, MP_2, MP_8, MP_9, desc(MP_3), desc(MP_4), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    tibble::column_to_rownames(var = "sample_id")
  mp_score_sorted <- mp_score_sorted[, -which(colnames(mp_score_sorted) == "ajcc_stage"), drop = FALSE]

  annotation_col <- annotation_col[rownames(mp_score_sorted), , drop = FALSE]

  pdf(output_path, width = 12, height = 3.5)
  pheatmap::pheatmap(
    t(mp_score_sorted),
    clustering_method = "ward.D2",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = "Tumor cell MP scores across sc-samples"
  )
  dev.off()
}

# ---- Main Function ----

#' Cluster tumor samples by MP ratio and MP score
#'
#' @param seurat_path Path to assigned MP Seurat RDS file.
#' @param output_dir Output directory for results.
run_analysis <- function(seurat_path, output_dir) {
  message("Loading single-cell data from: ", seurat_path)
  LUAD_Epi_assigned_MP <- readRDS(seurat_path)

  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  LUAD_Epi_assigned_MP <- reorder_mp_labels(LUAD_Epi_assigned_MP, old_order)

  message("Computing MP proportions per sample")
  count_table <- compute_mp_proportions(LUAD_Epi_assigned_MP)

  message("Building wide-format proportion matrix")
  df_wide <- count_table %>%
    tidyr::pivot_wider(
      id_cols = sample,
      names_from = MP,
      values_from = proportion
    ) %>%
    tibble::column_to_rownames(var = "sample")

  message("Building sample annotations")
  annotation_col <- build_sample_annotation(LUAD_Epi_assigned_MP, rownames(df_wide))

  ann_colors <- list(
    stage = c("Early" = "#1F77B4", "Advance" = "#D62728"),
    ajcc_stage = c(
      "I" = "#4DBBD5FF",
      "II" = "pink",
      "III" = "#F39B7FFF",
      "IV" = "#E64B35FF"
    )
  )

  figure_dir <- file.path(output_dir, "results_figure")

  proportion_path <- file.path(figure_dir, "03_tumor_ratio_cluster.pdf")
  message("Saving proportion heatmap to: ", proportion_path)
  plot_proportion_heatmap(df_wide, annotation_col, ann_colors, proportion_path)

  message("Computing per-sample MP scores")
  mp_score <- LUAD_Epi_assigned_MP@meta.data[, c("Sample", paste0("MP_", 1:9))]
  mp_score <- mp_score %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(dplyr::across(MP_1:MP_9, mean, na.rm = TRUE)) %>%
    tibble::column_to_rownames(var = "Sample")
  mp_score <- mp_score[rownames(df_wide), old_order]
  colnames(mp_score) <- paste0("MP_", 1:9)
  mp_score <- mp_score[, c("MP_1", "MP_2", "MP_8", "MP_9", "MP_3", "MP_4", "MP_5", "MP_6", "MP_7")]

  score_path <- file.path(figure_dir, "03_tumor_MP_socres_cluster.pdf")
  message("Saving MP score heatmap to: ", score_path)
  plot_score_heatmap(mp_score, annotation_col, ann_colors, score_path)

  message("Tumor ratio clustering complete")
}

# ---- Entry Point ----
if (!interactive()) {
  run_analysis(
    seurat_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    output_dir = RESULT_6_DIR
  )
}
