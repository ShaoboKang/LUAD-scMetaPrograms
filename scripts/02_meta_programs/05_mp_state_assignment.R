# ============================================================================
# Script: 05_mp_state_assignment
# Module: 02_meta_programs
# Description:
#   Load consensus MPs, re-order/re-name them, and visualize MP signatures
#   across assigned tumor cell states (heatmaps, violin plots).
#
# Prerequisites:
#   - 02_mp_signature_identification.R must have been run.
#   - Result_2/results_data/02_Malig_final_MP_top50.RData exists.
#   - Result_2/results_data/04_LUAD_Epi_assigned_MP_final.rds exists.
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
library(ggrepel)
library(patchwork)
library(dplyr)

# ---- Helper Functions ----

#' Load MP gene lists.
#' @param rdata_path Path to RData file.
#' @return Named list of MP gene vectors.
load_mp_list <- function(rdata_path) {
  message("Loading MP list from: ", rdata_path)
  load(rdata_path)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_list <- MP_list[old_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))
  mp_list
}

#' Load assigned MP Seurat object and standardize MP names.
#' @param rds_path Path to RDS file.
#' @return Seurat object with standardized MP column.
load_assigned_mp_data <- function(rds_path) {
  message("Loading assigned MP data from: ", rds_path)
  luad_epi <- readRDS(rds_path)
  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))
  luad_epi
}

#' Plot MP signature heatmap.
#' @param seurat_obj Seurat object.
#' @param mp_list MP gene list.
#' @param output_path Output PDF path.
#' @param title Plot title.
#' @param width PDF width.
#' @param height PDF height.
#' @param genes_show Optional vector of gene names to display.
plot_mp_heatmap <- function(seurat_obj, mp_list, output_path, title = "MP tumor cell state",
                            width = 15, height = 8, genes_show = NULL) {
  col <- c(
    "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
  )
  p <- DoHeatmap(
    seurat_obj,
    features = unique(unlist(mp_list)),
    group.by = "MP",
    label = FALSE,
    group.colors = col
  ) +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

  if (!is.null(genes_show)) {
    p <- p + scale_y_discrete(labels = function(x) ifelse(x %in% genes_show, x, ""))
  }

  p <- p + theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + labs(title = title) + coord_cartesian(clip = "off")

  pdf(file = output_path, width = width, height = height)
  print(p)
  dev.off()
  message("Saved MP heatmap to: ", output_path)
}

#' Plot violin plots of MP scores across cell states.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
#' @param width Width in inches.
#' @param height Height in inches.
plot_mp_violins <- function(seurat_obj, output_path, width = 16, height = 14) {
  col <- c(
    "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
  )
  plot_list <- list()
  for (i in paste0("MP_", 1:9)) {
    plot_list[[i]] <- VlnPlot(seurat_obj, features = i, pt.size = 0, group.by = "MP", cols = col) +
      geom_boxplot(width = 0.2, col = "black", fill = "white") +
      NoLegend()
  }
  p_combined <- wrap_plots(plot_list, nrow = 3, ncol = 3) +
    plot_annotation(
      title = "MP Scores in nine tumor cell state",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  ggsave(filename = output_path, plot = p_combined, width = width, height = height, dpi = 300)
  message("Saved MP violin plots to: ", output_path)
}

#' Plot MP score pheatmap with annotations.
#' @param seurat_obj Seurat object.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_mp_pheatmap <- function(seurat_obj, output_path, width = 15, height = 8) {
  annotation_col <- data.frame(cellname = colnames(seurat_obj), MP = seurat_obj$MP)
  annotation_col <- annotation_col[order(annotation_col$MP), ]
  annotation_col <- data.frame(row.names = annotation_col$cellname, MP = annotation_col$MP)

  data <- seurat_obj@meta.data[, paste0("MP_", 1:9)]
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  data <- data[, old_order]
  colnames(data) <- paste0("MP_", 1:9)
  data <- t(data[rownames(annotation_col), ])

  ann_colors <- list(
    MP = c(
      MP_1 = "#DC0000FF",
      MP_2 = "#E64B35FF",
      MP_3 = "#00A087FF",
      MP_4 = "#4DBBD5FF",
      MP_5 = "#3C5488FF",
      MP_6 = "#8491B4FF",
      MP_7 = "#F39B7FFF",
      MP_8 = "#7E6148FF",
      MP_9 = "#91D1C2FF"
    )
  )

  annotation_col$MP <- factor(annotation_col$MP, levels = unique(annotation_col$MP))

  bk <- c(seq(-4, -0.1, by = 0.01), seq(0, 4, by = 0.01))
  pdf(file = output_path, width = width, height = height)
  pheatmap::pheatmap(
    data, annotation_col = annotation_col, scale = "none", annotation_colors = ann_colors,
    cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, cluster_cols = FALSE,
    fontsize_row = 10, fontsize_col = 10, use_raster = FALSE, border_color = "black",
    color = c(
      colorRampPalette(colors = c("#091A8C", "white"))(length(bk) / 2),
      colorRampPalette(colors = c("white", "#B30000"))(length(bk) / 2)
    ),
    legend_breaks = seq(-4, 4, 2), breaks = bk
  )
  dev.off()
  message("Saved MP pheatmap to: ", output_path)
}

# ---- Main Function ----

#' Run MP state assignment visualization.
#' @param rdata_path Path to MP RData file.
#' @param rds_path Path to assigned MP Seurat object.
#' @param output_fig_dir Directory for output figures.
#' @return Seurat object used for plotting.
run_analysis <- function(
    rdata_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    rds_path = file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_list <- load_mp_list(rdata_path)
  luad_epi <- load_assigned_mp_data(rds_path)
  message("MP distribution:\n")
  print(table(luad_epi$MP))

  luad_epi <- NormalizeData(luad_epi) %>% ScaleData(features = unique(unlist(mp_list)))

  plot_mp_heatmap(luad_epi, mp_list, file.path(output_fig_dir, "09_MP_tumor_subcluster_signature.pdf"),
                  title = "MP tumor cell state")

  genes_show <- c(
    "TOP2A", "MKI67", "UBE2C", "UBE2T", "GINS2", "GMNN",
    "CXCL1", "CXCL2", "NFKBIA", "CCNL1",
    "CD74", "HLA-DRB1", "HLA-DRA", "HLA-DMA", "HLA-DPA1",
    "EIF3E", "EIF3F", "EIF1", "SLC25A6", "ATP6V1F",
    "CYB5A", "UQCR10", "UQCRB", "SLC25A5",
    "FXYD3", "TMSB10", "S100A10", "S100A11",
    "LAMC2", "ITGA2", "VEGFA",
    "ANGPTL4", "TGFBI", "ERO1A", "VEGFA", "PLOD2"
  )

  p024_t <- subset(luad_epi, Sample == "p024_T")
  plot_mp_heatmap(p024_t, mp_list, file.path(output_fig_dir, "10_p024_T_MP_tumor_subcluster_signature.pdf"),
                  title = "p024_T", genes_show = genes_show)

  plot_mp_violins(luad_epi, file.path(output_fig_dir, "15_MP_scores_in_tumor_vlnplot.pdf"))
  plot_mp_pheatmap(luad_epi, file.path(output_fig_dir, "16_MP_scores_in_tumor_heatmap.pdf"))

  message("MP state assignment visualization complete.")
  luad_epi
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP state assignment: ", e$message)
      stop(e)
    }
  )
}
