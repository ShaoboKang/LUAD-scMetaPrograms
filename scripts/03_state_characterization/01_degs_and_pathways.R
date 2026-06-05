# ============================================================================
# Script: 01_degs_and_pathways
# Module: 03_state_characterization
# Description: Differential expression and Hallmark GSVA for tumor cell states.
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
library(GSVA)
library(dplyr)
library(ggplot2)

# ---- Helper Functions ----

#' Load Seurat object and remap MP labels to standard order.
#' @param input_path Path to the Seurat RDS file.
#' @return Seurat object with remapped MP identities.
load_and_remap_seurat <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path)
  }

  message("[INFO] Loading Seurat object from: ", input_path)
  seurat_obj <- tryCatch({
    readRDS(input_path)
  }, error = function(e) {
    stop("Failed to read RDS file: ", conditionMessage(e))
  })

  seurat_obj$MP <- as.character(seurat_obj$MP)

  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  seurat_obj$MP <- unname(mp_map[seurat_obj$MP])
  seurat_obj$MP <- factor(seurat_obj$MP, levels = paste0("MP_", 1:9))

  return(seurat_obj)
}

#' Check for cross-lineage contamination with VlnPlot and DotPlot.
#' @param seurat_obj Seurat object.
#' @return List containing vln_plot and dot_plot (invisibly).
check_contamination <- function(seurat_obj) {
  message("[INFO] Checking for cross-lineage contamination...")

  vln_plot <- VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA"))

  validity_markers <- list(
    "Epithelial" = c("EPCAM", "KRT8", "KRT18", "KRT19", "MSLN", "MUC1"),
    "Immune" = c("PTPRC", "LST1", "TYROBP", "FCER1G", "TREM1", "S100A8"),
    "Plasma" = c("IGHG1", "IGHA1"),
    "Mycaf" = c("COL1A1", "COL1A2", "DCN", "LUM"),
    "Endothelial" = c("PECAM1", "VWF", "KDR", "EMCN")
  )

  dot_plot <- DotPlot(object = seurat_obj, features = validity_markers, group.by = "MP") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(title = "Expression(%)", order = 3),
      color = guide_colorbar(title = "Average Expression")
    ) +
    scale_color_gradientn(
      values = seq(0, 1, 0.2),
      colours = c("#f7f6fb", "#532b83")
    )

  print(vln_plot)
  print(dot_plot)

  invisible(list(vln_plot = vln_plot, dot_plot = dot_plot))
}

#' Find all markers for each MP cluster.
#' @param seurat_obj Seurat object with MP identities.
#' @return Data frame of DEGs.
find_degs <- function(seurat_obj) {
  message("[INFO] Finding DEGs with FindAllMarkers...")
  Idents(seurat_obj) <- "MP"

  degs <- FindAllMarkers(
    seurat_obj,
    min.pct = 0.25,
    logfc.threshold = 1,
    only.pos = TRUE
  )

  degs <- degs[which(degs$p_val_adj < 0.05), ]
  degs$group <- "adjust Pvalue < 0.05"
  degs$cluster <- factor(degs$cluster, levels = paste0("MP_", 1:9))

  message("[INFO] Found ", nrow(degs), " significant DEGs.")
  return(degs)
}

#' Plot DEGs as a jitter plot with top 30 labels per cluster.
#' @param degs Data frame of DEGs.
#' @param output_path Output PDF path.
#' @return NULL (invisibly).
plot_degs <- function(degs, output_path) {
  message("[INFO] Plotting DEGs to: ", output_path)

  # Group by cluster and label top 30 genes by avg_log2FC descending
  degs <- degs %>%
    group_by(cluster) %>%
    mutate(
      label = ifelse(
        row_number(desc(avg_log2FC)) <= 30,
        gene,
        NA
      )
    ) %>%
    ungroup()

  df_bg <- degs %>%
    group_by(cluster) %>%
    summarize(
      max_log2FC = max(avg_log2FC),
      min_log2FC = min(avg_log2FC),
      .groups = "drop"
    )

  p <- ggplot() +
    # Grey background on positive y-axis
    geom_col(
      data = df_bg,
      mapping = aes(cluster, max_log2FC),
      fill = "grey85", width = 0.8, alpha = 0.5
    ) +
    # Grey background on negative y-axis
    geom_col(
      data = df_bg,
      mapping = aes(cluster, min_log2FC),
      fill = "grey85", width = 0.8, alpha = 0.5
    ) +
    # Jittered data points
    geom_jitter(
      data = degs,
      mapping = aes(x = cluster, y = avg_log2FC, color = group),
      size = 1.5, width = 0.4, alpha = 0.7
    ) +
    # Group blocks at y = 0.6
    geom_col(
      data = df_bg,
      mapping = aes(x = cluster, y = 0.6, fill = cluster),
      width = 0.8
    ) +
    # Group labels inside blocks
    geom_text(
      data = df_bg,
      mapping = aes(x = cluster, y = 0.4, label = cluster),
      size = 5, color = "black", fontface = "bold"
    ) +
    scale_color_manual(values = c("#e42313", "#0061d5")) +
    scale_fill_manual(values = c(
      "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
      "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
    )) +
    # Label top genes
    geom_text_repel(
      data = degs,
      mapping = aes(x = cluster, y = avg_log2FC, label = label),
      max.overlaps = 10000,
      size = 3,
      segment.color = "black",
      show.legend = FALSE
    ) +
    theme_classic() +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 13, color = "black", face = "bold"),
      axis.line.y = element_line(color = "black", size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      legend.direction = "vertical",
      legend.justification = c(1, 0),
      legend.text = element_text(size = 13)
    ) +
    labs(x = "group", y = "Log2FoldChange", fill = NULL, color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    coord_cartesian(ylim = c(0.5, 6))

  pdf(file = output_path, width = 10, height = 6)
  print(p)
  dev.off()
  message("[INFO] DEG plot saved.")

  invisible(NULL)
}

#' Run GSVA with Hallmark gene sets on average expression per MP.
#' @param seurat_obj Seurat object.
#' @param degs Data frame of DEGs.
#' @param hallmark_path Path to Hallmark RData file.
#' @param output_path Output PDF path for heatmap.
#' @return GSVA matrix (invisibly).
run_gsva_hallmark <- function(seurat_obj, degs, hallmark_path, output_path) {
  if (!file.exists(hallmark_path)) {
    stop("Hallmark file not found: ", hallmark_path)
  }

  message("[INFO] Running Hallmark GSVA...")
  load(hallmark_path)

  exp <- AverageExpression(seurat_obj, group.by = "MP", slot = "data")[[1]]
  exp <- exp[unique(degs$gene), ]

  param <- gsvaParam(
    exprData = as.matrix(exp),
    geneSets = m_hallmark
  )

  gsva_matrix <- gsva(param)

  annotation_col <- data.frame(
    row.names = paste0("MP_", 1:9),
    MP = paste0("MP_", 1:9)
  )

  pdf(file = output_path)
  pheatmap::pheatmap(
    gsva_matrix,
    annotation_col = annotation_col,
    scale = "row",
    color = colorRampPalette(c("#091A8C", "white", "#B30000"))(100),
    show_colnames = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    cluster_rows = TRUE
  )
  dev.off()
  message("[INFO] Hallmark heatmap saved to: ", output_path)

  invisible(gsva_matrix)
}

# ---- Main Analysis Functions ----

#' Main analysis function for DEGs and pathways.
#' @param input_path Path to input Seurat RDS file.
#' @param output_dir Directory for output files.
#' @param msigdb_dir Directory containing MSigDB gene sets.
#' @return List of output file paths.
run_analysis <- function(input_path, output_dir, msigdb_dir) {
  # Ensure output subdirectories exist
  if (!dir.exists(file.path(output_dir, "results_data"))) {
    dir.create(file.path(output_dir, "results_data"), recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(file.path(output_dir, "results_figure"))) {
    dir.create(file.path(output_dir, "results_figure"), recursive = TRUE, showWarnings = FALSE)
  }

  # Step 1: Load and remap data
  seurat_obj <- load_and_remap_seurat(input_path)

  # Step 2: Contamination check plots
  check_contamination(seurat_obj)

  # Step 3: Find DEGs
  degs <- find_degs(seurat_obj)

  # Create DEG list (matching original behavior)
  deg_list <- split(degs$gene, degs$cluster)
  # saveRDS(deg_list, file = file.path(output_dir, "results_data", "01_tumor_DEG_list.rds"))

  # Step 4: Plot DEGs
  deg_plot_path <- file.path(output_dir, "results_figure", "01_tumor_state_DEGs.pdf")
  plot_degs(degs, deg_plot_path)

  # Step 5: Hallmark GSVA
  hallmark_path <- file.path(msigdb_dir, "m_hallmark.RData")
  gsva_plot_path <- file.path(output_dir, "results_figure", "samples_m_hallmark.pdf")
  run_gsva_hallmark(seurat_obj, degs, hallmark_path, gsva_plot_path)

  message("[INFO] Analysis complete.")
  return(list(
    deg_plot = deg_plot_path,
    gsva_heatmap = gsva_plot_path
  ))
}

# ---- Entry Point ----

if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_gsva_hallmark()
} else {
  message("[INFO] Script loaded. Call run_gsva_hallmark() to execute.")
}
