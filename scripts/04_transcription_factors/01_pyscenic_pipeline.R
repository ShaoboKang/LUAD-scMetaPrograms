#' pySCENIC Post-Processing Pipeline for LUAD MetaPrograms
# Script: 01_pyscenic_pipeline
#'
#' @description
#' Loads pySCENIC loom output, computes Regulon Specificity Scores (RSS),
#' visualizes AUC heatmaps, and correlates TF activity with MP signatures.
#'
#' @return Invisibly returns NULL. Writes figures and RDS to RESULT_4_DIR.
#'
#' @details
#' Run this script after completing the pySCENIC command-line workflow
#' (GRN, ctx, aucell) to generate the out_SCENIC.loom file.
#'
#' @examples
#' \dontrun{
#' run_pyscenic_pipeline()
#' }
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(ComplexHeatmap)
  library(gridExtra)
  library(viridis)
  library(ggrepel)
  library(cowplot)
})

#' Load and Reorder MP Labels
#'
#' @param rds_path Path to Seurat RDS file.
#' @return Seurat object with MP factor reordered.
load_mp_data <- function(rds_path) {
  obj <- readRDS(rds_path)
  obj$MP <- as.character(obj$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  obj$MP <- unname(mp_map[obj$MP])
  obj$MP <- factor(obj$MP, levels = paste0("MP_", 1:9))
  return(obj)
}

#' Load Regulon Data from Loom
#'
#' @param loom_path Path to SCENIC loom file.
#' @param seurat_obj Seurat object for cell name matching.
#' @return List with regulons, regulonAUC, regulonAucThresholds, sub_regulonAUC.
load_regulon_data <- function(loom_path, seurat_obj) {
  loom <- open_loom(loom_path)
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  regulonAucThresholds <- get_regulon_thresholds(loom)
  sub_regulonAUC <- regulonAUC[, match(colnames(seurat_obj), colnames(regulonAUC))]
  stopifnot(identical(colnames(sub_regulonAUC), colnames(seurat_obj)))
  close_loom(loom)
  list(
    regulons = regulons,
    regulonAUC = regulonAUC,
    regulonAucThresholds = regulonAucThresholds,
    sub_regulonAUC = sub_regulonAUC
  )
}

#' Plot RSS Heatmap
#'
#' @param rss RSS matrix.
#' @param out_pdf Output PDF path.
#' @return ggplot object.
plot_rss_heatmap <- function(rss, out_pdf) {
  rss <- rss[, paste0("MP_", 1:9)]
  rssPlot <- plotRSS(
    rss,
    labelsToDiscard = NULL,
    zThreshold = 1.5,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr = 0.01,
    varName = "cellType",
    col.low = "grey",
    col.mid = "#330066",
    col.high = "#330066",
    revCol = FALSE,
    verbose = TRUE
  )
  pdf(file = out_pdf, width = 4)
  print(rssPlot$plot)
  dev.off()
  return(rssPlot)
}

#' Plot TF RSS Rank per MP
#'
#' @param rss_df Data frame of RSS values.
#' @param celltypes Vector of cell type names.
#' @param out_pdf Output PDF path.
plot_tf_rss_rank <- function(rss_df, celltypes, out_pdf) {
  rssRanklist <- list()
  for (i in seq_along(celltypes)) {
    data_rank_plot <- data.frame(
      TF = rownames(rss_df),
      celltype = rss_df[, celltypes[i]],
      stringsAsFactors = FALSE
    )
    data_rank_plot <- na.omit(data_rank_plot)
    data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype, decreasing = TRUE), ]
    data_rank_plot$rank <- seq_len(nrow(data_rank_plot))

    p <- ggplot(data_rank_plot, aes(x = rank, y = celltype)) +
      geom_point(size = 3, shape = 16, color = "#1F77B4", alpha = 0.4) +
      geom_point(data = data_rank_plot[1:6, ], size = 3, color = "#DC050C") +
      theme_bw() +
      theme(
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(x = "Regulons Rank", y = "Specificity Score", title = celltypes[i]) +
      geom_text_repel(
        data = data_rank_plot[1:6, ],
        aes(label = TF),
        color = "black",
        size = 3,
        fontface = "italic",
        arrow = arrow(ends = "first", length = unit(0.01, "npc")),
        box.padding = 0.2,
        point.padding = 0.3,
        segment.color = "black",
        segment.size = 0.3,
        force = 1,
        max.iter = 3e3
      )
    rssRanklist[[i]] <- p
  }

  pdf(file = out_pdf, width = 24, height = 10)
  print(plot_grid(plotlist = rssRanklist, nrow = 2, ncol = 5, align = "vh"))
  dev.off()
}

#' Plot AUC Heatmap for Top TFs
#'
#' @param auc AUC matrix.
#' @param a_top5 Top 5 TF data frame.
#' @param annotation_col Column annotation data frame.
#' @param ann_colors Annotation colors list.
#' @param out_pdf Output PDF path.
plot_auc_heatmap <- function(auc, a_top5, annotation_col, ann_colors, out_pdf) {
  auc <- auc[rev(levels(a_top5$Topic)), rownames(annotation_col)]
  bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
  pdf(file = out_pdf, width = 15, height = 8)
  pheatmap::pheatmap(
    auc,
    annotation_col = annotation_col,
    scale = "row",
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    use_raster = FALSE,
    border_color = "black",
    color = c(
      colorRampPalette(colors = c("#091A8C", "white"))(length(bk) / 2),
      colorRampPalette(colors = c("white", "#B30000"))(length(bk) / 2)
    ),
    legend_breaks = seq(-4, 4, 2),
    breaks = bk
  )
  dev.off()
}

#' Plot Average Regulon Activity Heatmap
#'
#' @param regulonAUC Regulon AUC object.
#' @param annotation_col Column annotation data frame.
#' @param a_top5 Top 5 TF data frame.
#' @param out_pdf Output PDF path.
plot_avg_regulon_activity <- function(regulonAUC, annotation_col, a_top5, out_pdf) {
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
  regulonActivity_byCellType <- sapply(
    split(rownames(annotation_col), annotation_col$MP),
    function(cells) rowMeans(getAUC(regulonAUC)[, cells])
  )
  regulonActivity_byCellType <- regulonActivity_byCellType[rev(levels(a_top5$Topic)), ]
  regulonActivity_byCellType_Scaled <- scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE)

  pdf(file = out_pdf, width = 12, height = 4)
  Heatmap(
    regulonActivity_byCellType_Scaled,
    name = "Regulon activity",
    col = colorRamp2(c(-2, 0, 2), c("#091A8C", "white", "#B30000")),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  dev.off()
}

#' Plot Binary Regulon Activity Heatmap
#'
#' @param auc AUC matrix.
#' @param regulonAucThresholds Thresholds data frame.
#' @param annotation_col Column annotation data frame.
#' @param ann_colors Annotation colors list.
#' @param out_pdf Output PDF path.
plot_binary_regulon_activity <- function(auc, regulonAucThresholds, annotation_col, ann_colors, out_pdf) {
  regulonAucThresholds <- as.data.frame(regulonAucThresholds)
  regulonAucThresholds$Thresholds <- as.numeric(rownames(regulonAucThresholds))
  rownames(regulonAucThresholds) <- regulonAucThresholds$regulonAucThresholds

  regulonActivity <- auc
  for (i in rownames(regulonActivity)) {
    thresh <- regulonAucThresholds[i, "Thresholds"]
    regulonActivity[i, ] <- ifelse(regulonActivity[i, ] >= thresh, 1, 0)
  }

  pdf(file = out_pdf, width = 15, height = 8)
  pheatmap::pheatmap(
    regulonActivity,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    use_raster = FALSE,
    border_color = "black",
    color = c("white", "black")
  )
  dev.off()
}

#' Plot Correlation Between TF Activity and MP Signature
#'
#' @param auc AUC matrix.
#' @param mp_score MP signature matrix.
#' @param out_pdf Output PDF path.
plot_mp_auc_correlation <- function(auc, mp_score, out_pdf) {
  MP_auc_score <- rbind(mp_score, auc)
  r <- cor(t(MP_auc_score), method = "pearson", use = "pairwise.complete.obs")
  r <- r[1:9, -c(1:9)]

  pdf(file = out_pdf, width = 12, height = 4)
  Heatmap(
    r,
    name = "Correlation",
    col = colorRamp2(c(-0.5, 0, 0.5), c("#091A8C", "white", "#B30000")),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  dev.off()
}

#' Run Full pySCENIC Post-Processing Pipeline
#'
#' @return Invisibly returns NULL.
run_pyscenic_pipeline <- function() {
  message("[01_pyscenic_pipeline] Starting pipeline...")

  mp_rds <- file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds")
  loom_file <- file.path(RESULT_4_DIR, "scenic", "output", "out_SCENIC.loom")

  LUAD_Epi_assigned_MP <- load_mp_data(mp_rds)

  regulon_data <- load_regulon_data(loom_file, LUAD_Epi_assigned_MP)
  regulons <- regulon_data$regulons
  regulonAUC <- regulon_data$regulonAUC
  regulonAucThresholds <- regulon_data$regulonAucThresholds
  sub_regulonAUC <- regulon_data$sub_regulonAUC

  LUAD_Epi_assigned_MP@meta.data <- cbind(
    LUAD_Epi_assigned_MP@meta.data,
    t(sub_regulonAUC@assays@data$AUC)
  )

  rss <- calcRSS(
    AUC = getAUC(sub_regulonAUC),
    cellAnnotation = LUAD_Epi_assigned_MP$MP
  )

  rssPlot <- plot_rss_heatmap(
    rss,
    file.path(RESULT_4_DIR, "results_figure", "06_special_TF.pdf")
  )

  B_rss <- as.data.frame(rss)
  celltype <- paste0("MP_", 1:9)
  plot_tf_rss_rank(
    B_rss,
    celltype,
    file.path(RESULT_4_DIR, "results_figure", "01_TF_rss.pdf")
  )

  rssPlot$df$cellType <- factor(rssPlot$df$cellType, levels = paste0("MP_", 1:9))
  a_top5 <- rssPlot$df %>%
    group_by(cellType) %>%
    arrange(desc(RSS), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()

  auc <- getAUC(sub_regulonAUC)
  annotation_col <- data.frame(
    cellname = colnames(LUAD_Epi_assigned_MP),
    MP = LUAD_Epi_assigned_MP$MP
  )
  annotation_col <- annotation_col[order(annotation_col$MP), ]
  annotation_col <- data.frame(row.names = annotation_col$cellname, MP = annotation_col$MP)
  annotation_col$MP <- factor(annotation_col$MP, levels = unique(annotation_col$MP))

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

  plot_auc_heatmap(
    auc,
    a_top5,
    annotation_col,
    ann_colors,
    file.path(RESULT_4_DIR, "results_figure", "02_Auc_scores_in_tumor_heatmap.pdf")
  )

  plot_avg_regulon_activity(
    regulonAUC,
    annotation_col,
    a_top5,
    file.path(RESULT_4_DIR, "results_figure", "03_Auc_scores_avg_in_tumor_heatmap.pdf")
  )

  plot_binary_regulon_activity(
    auc,
    regulonAucThresholds,
    annotation_col,
    ann_colors,
    file.path(RESULT_4_DIR, "results_figure", "04_Auc_binary_in_tumor_heatmap.pdf")
  )

  MP_score <- t(LUAD_Epi_assigned_MP@meta.data[colnames(auc), paste0("MP_", 1:9)])
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  MP_score <- MP_score[old_order, ]
  rownames(MP_score) <- paste0("MP_", 1:9)
  stopifnot(all(colnames(auc) == colnames(MP_score)))

  plot_mp_auc_correlation(
    auc,
    MP_score,
    file.path(RESULT_4_DIR, "results_figure", "05_Auc_MP_cor_in_tumor_heatmap.pdf")
  )

  saveRDS(a_top5, file = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"))

  message("[01_pyscenic_pipeline] Completed successfully.")
  invisible(NULL)
}

if (!interactive()) {
  run_pyscenic_pipeline()
}
