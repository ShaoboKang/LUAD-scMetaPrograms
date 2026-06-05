# ============================================================================
# Script: 05_clustergvis_all_mp_cells
# Module: 07_trajectory
# Description:
#   ClusterGVis visualization of TF activities along Monocle3 pseudotime.
#   Produces both line-and-heatmap plot and annotated heatmap with
#   branch/pseudotime color bars. Matches original ClusterGVis.R logic.
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

#' Run ClusterGVis analysis for all MP cells
#'
#' @param mp_list_path Path to MP top50 RData.
#' @param tf_path Path to top5 TF RDS.
#' @param tf_pseudotime_path Path to TF_activate_pseudotime RDS.
#' @param output_fig_dir Directory for figures.
#' @return List with cluster objects and annotation.
run_clustergvis_all_mp_cells <- function(
    mp_list_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    tf_path = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"),
    tf_pseudotime_path = file.path(RESULT_7_DIR, "results_data", "TF_activate_pseudotime.rds"),
    output_fig_dir = file.path(RESULT_7_DIR, "results_figure")
) {
  library(org.Hs.eg.db)
  library(ClusterGVis)
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(ComplexHeatmap)

  if (!dir.exists(output_fig_dir)) {
    dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
  }

  message("Loading data...")
  load(mp_list_path)

  tf_names <- readRDS(tf_path)
  tf_names <- unique(tf_names$Topic)

  tf_activate_pseudotime <- readRDS(tf_pseudotime_path)
  tf_activate_pseudotime <- tf_activate_pseudotime[
    !tf_activate_pseudotime$branch %in% c("7", "5"), ]
  tf_activate_pseudotime <- tf_activate_pseudotime[
    order(tf_activate_pseudotime$pseudotime), ]
  tf_activate_pseudotime$branch <- factor(
    tf_activate_pseudotime$branch, levels = c("1", "2", "3", "4", "6")
  )

  lineage_1 <- tf_activate_pseudotime[tf_activate_pseudotime$branch %in% c("2", "1", "4"), ]
  lineage_1 <- lineage_1[order(lineage_1$pseudotime, decreasing = TRUE), ]
  lineage_2 <- tf_activate_pseudotime[tf_activate_pseudotime$branch %in% c("6", "3"), ]
  lineage_2 <- lineage_2[order(lineage_2$pseudotime), ]
  lineage <- rbind(lineage_1, lineage_2)

  mat <- as.data.frame(t(lineage[, 4:140]))
  mat <- t(apply(mat, 1, function(x) {
    stats::smooth.spline(x, df = 3)$y
  }))
  mat <- t(apply(mat, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  colnames(mat) <- lineage$cell_id

  set.seed(123)
  ck <- clusterData(as.data.frame(mat), seed = 123,
                    cluster.method = "kmeans", cluster.num = 3)

  pdf(file.path(output_fig_dir, "monocle3_builtin_anno.pdf"),
      height = 12, width = 10, onefile = FALSE)
  visCluster(
    object = ck, plot.type = "both",
    add.sampleanno = FALSE, line.side = "left",
    cluster.order = c(3, 1, 2),
    markGenes = tf_names, cluster_columns = FALSE
  )
  dev.off()

  branch_colors <- c(
    "1" = "#E64B35FF", "2" = "#4DBBD5FF", "3" = "#91D1C2FF",
    "4" = "#00A087FF", "6" = "#F39B7FFF"
  )
  pseudotime_col_fun <- colorRamp2(
    seq(min(lineage$pseudotime), max(lineage$pseudotime), length = 100),
    viridisLite::viridis(100)
  )

  ha <- HeatmapAnnotation(
    Cluster = lineage$branch,
    Pseudotime = lineage$pseudotime,
    col = list(Cluster = branch_colors, Pseudotime = pseudotime_col_fun),
    annotation_name_side = "left"
  )

  rownames(mat) <- gsub("\\(.*\\)", "", rownames(mat))
  ck1 <- clusterData(as.data.frame(mat), seed = 123,
                     cluster.method = "kmeans", cluster.num = 3)

  anno <- enrichCluster(object = ck1,
                        OrgDb = org.Hs.eg.db,
                        type = "BP", fromType = "SYMBOL", organism = "hsa",
                        pvalueCutoff = 0.05, topn = 100)

  pdf(file.path(output_fig_dir, "monocle3_with_heatmapAnnotation.pdf"),
      height = 12, width = 10, onefile = FALSE)
  visCluster(
    object = ck, plot.type = "heatmap",
    HeatmapAnnotation = ha,
    sample.order = lineage$cell_id,
    line.side = "left", cluster.order = c(2, 1, 3),
    markGenes = tf_names, cluster_columns = FALSE
  )
  dev.off()

  message("ClusterGVis all MP cells complete.")
  return(list(ck = ck, ck1 = ck1, anno = anno))
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_clustergvis_all_mp_cells()
} else {
  message("[INFO] Script loaded. Call run_clustergvis_all_mp_cells() to execute.")
}
