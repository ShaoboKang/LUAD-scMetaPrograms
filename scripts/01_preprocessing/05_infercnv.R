# ============================================================================
# Script: 05_infercnv
# Module: 01_preprocessing
# Description:
#   Prepare inputs for inferCNV, run copy-number analysis on epithelial cells,
#   compute CNV scores, and visualize results.
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
  library(infercnv)
  library(copykat)
  library(ggplot2)

# ---- Helper Functions ----

#' Load final epithelial annotation.
#' @param path Path to .rds file.
#' @return Seurat object.
load_epithelial_data <- function(path) {
  message("Loading epithelial data from: ", path)
  readRDS(path)
}

#' Build inferCNV labels: normal epithelial clusters from normal tissue.
#' @param seurat_obj Seurat object.
#' @param normal_clusters Vector of cluster IDs considered normal.
#' @return Seurat object with infercnv_label column.
build_infercnv_labels <- function(seurat_obj, normal_clusters = c(1, 6, 7, 10, 22)) {
  seurat_obj$infercnv_label <- paste0(seurat_obj$seurat_clusters, "_T")
  for (cluster in normal_clusters) {
    idx <- which(seurat_obj$seurat_clusters == cluster & seurat_obj$Tissue == "Normal")
    seurat_obj$infercnv_label[idx] <- paste0(seurat_obj$seurat_clusters[idx], "_N")
  }
  seurat_obj
}

#' Load gene position information.
#' @param path Path to gene position file.
#' @return Data frame with gene positions.
load_gene_info <- function(path) {
  gene_info <- read.table(path)
  colnames(gene_info) <- c("gene_name", "chr", "start", "end")
  gene_info <- gene_info[!duplicated(gene_info[, 1]), ]
  allowed_chrs <- paste0("chr", 1:22)
  gene_info <- gene_info[gene_info$chr %in% allowed_chrs, ]
  gene_info$chr <- as.numeric(sub("^chr", "", gene_info$chr))
  gene_info <- gene_info[with(gene_info, order(chr, start)), ]
  gene_info
}

#' Prepare and write inferCNV input files.
#' @param expr_matrix Expression matrix (genes x cells).
#' @param gene_info Gene position data frame.
#' @param group_info Group annotation data frame.
#' @param output_dir Directory for output files.
#' @return List of written file paths.
write_infercnv_inputs <- function(expr_matrix, gene_info, group_info, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  exp_file <- file.path(output_dir, "expFile.txt")
  group_file <- file.path(output_dir, "groupFiles.txt")
  gene_file <- file.path(output_dir, "geneFile.txt")
  write.table(expr_matrix, file = exp_file, sep = "\t", quote = FALSE)
  write.table(group_info, file = group_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(gene_info, file = gene_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  list(exp = exp_file, group = group_file, gene = gene_file)
}

#' Run inferCNV analysis.
#' @param exp_file Expression file path.
#' @param group_file Group file path.
#' @param gene_file Gene order file path.
#' @param ref_groups Character vector of reference group names.
#' @param out_dir Output directory for inferCNV.
#' @return inferCNV object.
run_infercnv <- function(exp_file, group_file, gene_file, ref_groups, out_dir) {
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = exp_file,
    annotations_file = group_file,
    delim = "\t",
    gene_order_file = gene_file,
    ref_group_names = ref_groups
  )
  infercnv_obj2 <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    hclust_method = "ward.D2",
    plot_steps = FALSE,
    HMM = FALSE,
    denoise = TRUE,
    scale_data = TRUE
  )
  infercnv_obj2
}

#' Compute CNV scores from inferCNV expression data.
#' @param infercnv_obj inferCNV result object.
#' @param group_info Group annotation data frame.
#' @return Data frame with CNV scores.
compute_cnv_scores <- function(infercnv_obj, group_info) {
  expr <- infercnv_obj@expr.data
  expr2 <- (expr - 1)^2
  cnv_score <- as.data.frame(colMeans(expr2))
  colnames(cnv_score) <- "CNV_score"
  cnv_score$cell <- rownames(cnv_score)
  cnv_score$cluster <- group_info$cluster[match(cnv_score$cell, group_info$cellname)]
  cnv_score$CNV_score[which(cnv_score$CNV_score > 0.02)] <- 0.02
  cnv_score
}

#' Plot CNV scores by cluster.
#' @param cnv_score Data frame with CNV_score, cluster, and tissue columns.
#' @param output_path Output PDF path.
plot_cnv_scores <- function(cnv_score, output_path) {
  p <- ggplot(cnv_score, aes(cluster, CNV_score, fill = tissue)) +
    geom_violin(alpha = 0.4) +
    stat_boxplot(
      geom = "errorbar",
      position = position_dodge(width = 0.1),
      width = 0.1
    ) +
    geom_boxplot(alpha = 0.5, outlier.size = 0) +
    scale_fill_manual(
      values = c(
        "Normal_epi_from_normal" = "#fce38a",
        "Normal_epi_from_tumor" = "#1F77B4",
        "Tumor_epi_from_tumor" = "#D62728"
      )
    ) +
    theme_bw() +
    labs(x = "Epi cluster", y = "CNV Score") +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  pdf(file = output_path, width = 10, height = 6)
  print(p)
  dev.off()
  message("Saved CNV score plot to: ", output_path)
}

# ---- Main Function ----

#' Run inferCNV analysis pipeline.
#' @param input_path Path to final epithelial Seurat object.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return List containing inferCNV object and CNV scores.
run_analysis <- function(
    input_path = file.path(RESULT_1_DIR, "results_data", "08_LUAD_epi_anno_final.rds"),
    output_data_dir = file.path(RESULT_1_DIR, "results_data", "infercnv"),
    output_fig_dir = file.path(RESULT_1_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  luad_epi <- load_epithelial_data(input_path)
  luad_epi <- build_infercnv_labels(luad_epi, normal_clusters = c(1, 6, 7, 10, 22))

  dat <- GetAssayData(luad_epi, layer = "counts")
  dat <- dat[which(rowSums(dat) > 0), ]
  colnames(dat) <- gsub("-", "_", colnames(dat))
  colnames(dat) <- gsub("_", "", colnames(dat))

  gene_info <- load_gene_info(file.path(output_data_dir, "gencode_v32_gene_pos_gene_name.txt"))
  gene_info <- gene_info[match(rownames(dat), gene_info$gene_name), ]
  gene_info <- na.omit(gene_info)

  group_info <- data.frame(
    cellname = colnames(dat),
    cluster = luad_epi$infercnv_label,
    tissue = luad_epi$Tissue
  )

  dat <- dat[gene_info$gene_name, ]
  message("Expression matrix dimensions: ", paste(dim(dat), collapse = " x "))

  inputs <- write_infercnv_inputs(dat, gene_info, group_info, output_data_dir)

  infercnv_obj2 <- run_infercnv(
    inputs$exp, inputs$group, inputs$gene,
    ref_groups = c("1_N", "6_N", "7_N", "10_N", "22_N"),
    out_dir = file.path(output_data_dir, "infercnv_output")
  )

  # Reorder observation groups
  current_names <- names(infercnv_obj2@observation_grouped_cell_indices)
  last_groups <- rev(c("1_T", "6_T", "7_T", "10_T", "22_T"))
  other_groups <- setdiff(current_names, last_groups)
  new_order <- c(other_groups, last_groups)
  infercnv_obj2@observation_grouped_cell_indices <- infercnv_obj2@observation_grouped_cell_indices[new_order]

  plot_cnv(infercnv_obj2, out_dir = output_data_dir, output_format = "pdf")
  saveRDS(infercnv_obj2, file = file.path(output_data_dir, "infercnv_obj2.rds"))

  cnv_score <- compute_cnv_scores(infercnv_obj2, group_info)

  normal_epi_from_normal <- c("1_N", "6_N", "7_N", "10_N", "22_N")
  normal_epi_from_tumor <- c("1_T", "6_T", "7_T", "10_T", "22_T")
  cnv_score$tissue <- ifelse(
    cnv_score$cluster %in% normal_epi_from_normal, "Normal_epi_from_normal",
    ifelse(cnv_score$cluster %in% normal_epi_from_tumor, "Normal_epi_from_tumor", "Tumor_epi_from_tumor")
  )

  # Sort by tissue and median CNV score
  cluster_median <- cnv_score %>%
    group_by(tissue, cluster) %>%
    summarise(CNV_median = median(CNV_score), .groups = "drop")

  cnv_score <- cnv_score %>%
    left_join(cluster_median, by = c("tissue", "cluster")) %>%
    arrange(tissue, CNV_median) %>%
    select(-CNV_median)

  cnv_score$cluster <- factor(cnv_score$cluster, levels = unique(cnv_score$cluster))

  plot_cnv_scores(cnv_score, file.path(output_fig_dir, "15_cnvscores_epi_cluster.pdf"))

  luad_epi$CNV_score <- cnv_score$CNV_score[match(colnames(luad_epi), cnv_score$cell)]

  pdf(file = file.path(output_fig_dir, "16_umap_cnvscores.pdf"), width = 8, height = 6)
  print(FeaturePlot(luad_epi, reduction = "umap.harmony", features = "CNV_score") +
    scale_color_viridis_c(option = "D"))
  dev.off()
  message("Saved CNV UMAP to: ", file.path(output_fig_dir, "16_umap_cnvscores.pdf"))

  saveRDS(cnv_score, file = file.path(output_data_dir, "09_CNV_score.rds"))
  message("inferCNV analysis complete.")

  invisible(list(infercnv_obj = infercnv_obj2, cnv_score = cnv_score))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in inferCNV analysis: ", e$message)
      stop(e)
    }
  )
}
