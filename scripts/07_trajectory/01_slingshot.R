# ============================================================================
# Script: 01_slingshot
# Module: 07_trajectory
# Description: Slingshot trajectory inference and pseudotime visualization.
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

#' Plot pseudotime density by cluster or MP
#'
#' @param seurat_obj A Seurat object with pseudotime metadata.
#' @param lineage Column name for lineage pseudotime.
#' @param cluster_label Metadata column for grouping.
#' @param color Color palette vector.
#' @return A ggplot object.
plot_pseudotime_density <- function(seurat_obj, lineage, cluster_label, color) {
  df <- data.frame(seurat_obj[[cluster_label]], seurat_obj[[lineage]])
  colnames(df) <- c("celltype", "lineage")
  df <- na.omit(df)
  p <- ggplot(df, aes(x = lineage, fill = celltype)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = color)
  return(p)
}

# ---- Main Function ----

#' Run Slingshot trajectory analysis
#'
#' Integrates Monocle3 clustering with Slingshot to infer trajectories,
#' plots UMAP with lineages, pseudotime heatmaps, and cell density along
#' each lineage.
#'
#' @param cds_path Path to Monocle3 CDS RDS.
#' @param sc_path Path to Seurat RDS.
#' @param output_dir Directory for outputs.
#' @return Slingshot SingleCellExperiment object.
run_analysis <- function(
    cds_path = file.path(RESULT_7_DIR, "results_data", "03_cds_monocle3.rds"),
    sc_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds"),
    output_dir = file.path(RESULT_7_DIR, "results_data")
) {
  library(slingshot)
  library(dplyr)
  library(Seurat)
  library(scales)
  library(fields)
  library(RColorBrewer)
  library(dittoSeq)
  library(ggplot2)

  message("Loading Monocle3 CDS and Seurat object...")
  cds <- readRDS(cds_path)
  luad_epi <- readRDS(sc_path)

  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))

  luad_epi$cluster <- paste0("cluster_", cds@clusters@listData[["UMAP"]][["clusters"]])
  luad_epi <- subset(luad_epi, cluster != "cluster_7")

  # Plot UMAP by cluster
  p_umap <- DimPlot(luad_epi, reduction = "umap.harmony", group.by = "cluster", label = TRUE)
  print(p_umap)

  sce <- as.SingleCellExperiment(luad_epi, assay = "RNA")

  message("Running Slingshot...")
  sce_slingshot <- slingshot(
    sce,
    reducedDim = "UMAP.HARMONY",
    clusterLabels = "cluster",
    start.clus = "cluster_2",
    end.clus = NULL,
    extend = "n"
  )

  SlingshotDataSet(sce_slingshot)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sce_slingshot, file = file.path(output_dir, "02_sce_slingshot1.rds"))

  # Plot cell trajectory on UMAP
  fig_dir <- file.path(RESULT_7_DIR, "results_figure")
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  sce_slingshot$MP <- factor(sce_slingshot$MP, levels = paste0("MP_", 1:9))
  sce_slingshot$cluster <- factor(sce_slingshot$cluster, levels = paste0("cluster_", 1:6))

  celltypes <- levels(sce_slingshot$cluster)
  col <- c("#E64B35FF", "#4DBBD5FF", "#91D1C2FF", "#00A087FF",
           "#3C5488FF", "#F39B7FFF")
  pal <- setNames(col[seq_along(celltypes)], celltypes)

  pdf(file.path(fig_dir, "05_slingshot.pdf"), width = 6, height = 6)
  plot(reducedDims(sce_slingshot)$UMAP.HARMONY,
       col = pal[sce_slingshot$cluster], cex = 0.3, pch = 16, asp = 1)
  lines(SlingshotDataSet(sce_slingshot), lwd = 2, col = brewer.pal(3, "Set1"))
  legend("topright", legend = names(pal), col = pal, pch = 16, cex = 0.6, box.lwd = 0.5)
  legend("bottomright",
         legend = paste0("lineage", 1:3),
         col = unique(brewer.pal(3, "Set1")), pch = 16, cex = 0.6, box.lwd = 0.5)
  title("Slingshot Trajectory Inference", cex.main = 1.2, font.main = 2)
  dev.off()

  # Pseudotime plots for 3 lineages
  pdf(file.path(fig_dir, "06_slingshot_pesudotime.pdf"), width = 18, height = 6)
  par(mfrow = c(1, 3), mar = c(5, 4, 4, 1) + 0.1)

  color_func <- function(pt) {
    rev(colorRampPalette(brewer.pal(11, "Spectral")[-6])(100))[
      cut(pt, breaks = 100, include.lowest = TRUE, na.rm = TRUE)]
  }

  for (i in 1:3) {
    pt <- sce_slingshot[[paste0("slingPseudotime_", i)]]
    plot(reducedDims(sce_slingshot)$UMAP.HARMONY,
         col = ifelse(is.na(pt), "gray", color_func(pt)),
         pch = 16, cex = 0.5,
         main = paste("Lineage", i),
         xlab = "UMAP1", ylab = "UMAP2")
    lines(SlingshotDataSet(sce_slingshot), lwd = 2, col = brewer.pal(3, "Set1"))
    legend("topleft",
           legend = paste0("Lineage ", i),
           col = brewer.pal(3, "Set1")[i],
           lty = 1, lwd = 2, cex = 0.8)
    if (i == 3) {
      image.plot(
        zlim = range(c(sce_slingshot$slingPseudotime_1,
                       sce_slingshot$slingPseudotime_2,
                       sce_slingshot$slingPseudotime_3), na.rm = TRUE),
        col = rev(colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)),
        legend.only = TRUE,
        horizontal = FALSE,
        legend.args = list(text = "Pseudotime", side = 4, line = 2.5, cex = 0.8),
        smallplot = c(0.85, 0.88, 0.65, 0.85),
        axis.args = list(cex.axis = 0.7)
      )
    }
  }
  par(mfrow = c(1, 1))
  dev.off()

  # Cell abundance density along pseudotime
  pseudotime <- slingPseudotime(sce_slingshot) %>% as.data.frame()
  lineages <- colnames(pseudotime)

  luad_epi <- AddMetaData(object = luad_epi,
                          metadata = pseudotime,
                          col.name = lineages)

  p1 <- plot_pseudotime_density(luad_epi, lineage = "Lineage1",
                                 cluster_label = "cluster", color = col) +
    theme(axis.title.x = element_blank())
  p2 <- plot_pseudotime_density(luad_epi, lineage = "Lineage2",
                                 cluster_label = "MP", color = col) +
    theme(axis.title.x = element_blank())

  print(p1)
  print(p2)

  message("Slingshot analysis complete.")
  return(sce_slingshot)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
