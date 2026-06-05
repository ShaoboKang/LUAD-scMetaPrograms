# ============================================================================
# Script: 02_bayesprism
# Module: 06_deconvolution
# Description: Run BayesPrism deconvolution using MP cell states.
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

#' Prepare single-cell count matrix from a Seurat object
#'
#' @param seurat_obj A Seurat object with counts assay.
#' @return Transposed count matrix.
prepare_sc_counts <- function(seurat_obj) {
  counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
  return(t(counts))
}

#' Create boxplots comparing MP fractions across clinical stages
#'
#' @param theta_df Data frame of cell type fractions with group column.
#' @return A patchwork object of boxplots.
plot_theta_boxplots <- function(theta_df) {
  target_columns <- paste0("MP_", 1:9)

  boxplot_list <- list()
  for (i in target_columns) {
    boxplot_list[[i]] <- ggpubr::ggboxplot(
      theta_df,
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

#' Run BayesPrism deconvolution analysis and save outputs
#'
#' @param seurat_path Path to assigned MP Seurat RDS file.
#' @param bulk_path Path to bulk clinical RData file.
#' @param output_dir Output directory for results.
run_analysis <- function(seurat_path, bulk_path, output_dir) {
  message("Loading bulk data from: ", bulk_path)
  load(bulk_path)
  bk_dat <- t(mRNA_exp)

  message("Loading single-cell data from: ", seurat_path)
  LUAD_Epi_assigned_MP <- readRDS(seurat_path)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  LUAD_Epi_assigned_MP <- reorder_mp_labels(LUAD_Epi_assigned_MP, old_order)

  message("Preparing single-cell count matrix")
  sc_dat <- prepare_sc_counts(LUAD_Epi_assigned_MP)

  cell_state_labels <- LUAD_Epi_assigned_MP$MP
  cell_type_labels <- LUAD_Epi_assigned_MP$MP

  rm(LUAD_Epi_assigned_MP)
  gc()

  data_dir <- file.path(output_dir, "results_data")
  figure_dir <- file.path(output_dir, "results_figure")

  cor_path <- file.path(data_dir, "celltype_cor_conbin.pdf")
  message("Plotting cell-state correlation to: ", cor_path)
  pdf(file = cor_path, width = 16, height = 16)
  BayesPrism::plot.cor.phi(
    input = sc_dat,
    input.labels = cell_state_labels,
    title = "cell state correlation",
    cexRow = 0.2,
    cexCol = 0.2,
    margins = c(2, 2)
  )
  dev.off()

  message("Cleaning up gene sets")
  sc_dat_filtered <- BayesPrism::cleanup.genes(
    input = sc_dat,
    input.type = "count.matrix",
    species = "hs",
    gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
    exp.cells = 5
  )

  message("Comparing bulk and single-cell expression")
  BayesPrism::plot.bulk.vs.sc(sc.input = sc_dat_filtered, bulk.input = bk_dat)

  message("Selecting protein-coding genes")
  sc_dat_filtered_pc <- BayesPrism::select.gene.type(sc_dat_filtered, gene.type = "protein_coding")

  message("Building BayesPrism model")
  my_prism <- BayesPrism::new.prism(
    reference = sc_dat_filtered_pc,
    mixture = bk_dat,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )

  message("Running BayesPrism (this may take a while)")
  bp_res <- BayesPrism::run.prism(prism = my_prism, n.cores = 32)

  bp_path <- file.path(data_dir, "BayesPrism.rds")
  message("Saving BayesPrism results to: ", bp_path)
  saveRDS(bp_res, file = bp_path)

  message("Extracting cell-type fractions")
  theta <- BayesPrism::get.fraction(
    bp = bp_res,
    which.theta = "final",
    state.or.type = "type"
  )
  theta <- as.data.frame(theta)
  theta$group <- clinical$stage
  theta$group <- factor(theta$group, levels = c("normal", "early", "advance"))

  message("Creating boxplots")
  boxplot_figure <- plot_theta_boxplots(theta)

  figure_path <- file.path(figure_dir, "02_BayesPrism.pdf")
  message("Saving figure to: ", figure_path)
  pdf(file = figure_path, width = 12, height = 12)
  print(boxplot_figure)
  dev.off()

  message("Plotting heatmap")
  pheatmap::pheatmap(
    t(cibersort),
    clustering_method = "ward.D2",
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = "Tumor cell state ratio across sc-samples"
  )

  message("BayesPrism analysis complete")
}

# ---- Entry Point ----
if (!interactive()) {
  run_analysis(
    seurat_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA", "TCGA_sur_data_stage.Rdata"),
    output_dir = RESULT_6_DIR
  )
}
