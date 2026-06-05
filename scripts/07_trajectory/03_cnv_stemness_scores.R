# ============================================================================
# Script: 03_cnv_stemness_scores
# Module: 07_trajectory
# Description: CNV score violin plots, CytoTRACE stemness analysis, and MP
#   proportion pie charts across clusters.
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

#' Compute median score per MP and reorder factor levels
#'
#' @param df Data frame with MP and Score columns.
#' @return Data frame with MP reordered by median score.
reorder_by_median <- function(df) {
  cluster_median <- df %>%
    group_by(MP) %>%
    summarise(score_median = median(Score, na.rm = TRUE), .groups = "drop")

  df <- df %>%
    left_join(cluster_median, by = "MP") %>%
    arrange(score_median) %>%
    select(-score_median)

  df$MP <- factor(df$MP, levels = unique(df$MP))
  return(df)
}

#' Plot violin and boxplot for scores by MP
#'
#' @param df Data frame with MP and Score columns.
#' @param y_label Y-axis label.
#' @param colors Color palette.
#' @return A ggplot object.
plot_score_violin <- function(df, y_label, colors) {
  ggplot(df, aes(x = MP, y = Score, fill = MP)) +
    geom_violin(alpha = 0.4) +
    stat_boxplot(
      geom = "errorbar",
      position = position_dodge(width = 0.1),
      width = 0.1
    ) +
    geom_boxplot(alpha = 0.5, outlier.size = 0, size = 0.3, width = 0.3) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    labs(x = "MP", y = y_label) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ---- Main Function ----

#' Run CNV and stemness score analysis
#'
#' Computes CNV violin plots, CytoTRACE stemness scores, and pie charts
#' showing MP proportions per cluster.
#'
#' @param sc_path Path to Seurat RDS.
#' @param cds_path Path to Monocle3 CDS RDS.
#' @param output_data_dir Directory for result RDS files.
#' @param output_fig_dir Directory for result figures.
#' @return List containing CytoTRACE results.
run_analysis <- function(
    sc_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds"),
    cds_path = file.path(RESULT_7_DIR, "results_data", "03_cds_monocle3.rds"),
    output_data_dir = file.path(RESULT_7_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_7_DIR, "results_figure")
) {
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggsci)
  library(dplyr)
  library(CytoTRACE)
  library(Seurat)

  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading Monocle3 CDS...")
  cds <- readRDS(cds_path)

  # CNV score data frame
  cnv_score <- data.frame(
    cellnames = colnames(cds),
    MP = cds@colData@listData[["MP"]],
    CNV_score = cds@colData@listData[["CNV_score"]]
  )

  cnv_score <- reorder_by_median(cnv_score %>% rename(Score = CNV_score))
  cnv_score <- cnv_score %>% rename(CNV_score = Score)

  cnv_colors <- c(
    "#3C5488FF", "#8491B4FF", "#4DBBD5FF", "#00A087FF",
    "#F39B7FFF", "#91D1C2FF", "#DC0000FF", "#E64B35FF", "#7E6148FF"
  )

  pdf(file.path(output_fig_dir, "CNV_score.pdf"), width = 8, height = 4)
  print(plot_score_violin(cnv_score %>% rename(Score = CNV_score), "CNV Score", cnv_colors))
  dev.off()

  # Stemness analysis with CytoTRACE
  message("Loading Seurat object for CytoTRACE...")
  luad_epi <- readRDS(sc_path)

  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))

  exp_mat <- as.matrix(GetAssayData(luad_epi, slot = "counts"))
  exp_mat <- exp_mat[apply(exp_mat > 0, 1, sum) >= 5, ]

  message("Running CytoTRACE...")
  cyto_results <- CytoTRACE(exp_mat, ncores = 10)
  saveRDS(cyto_results, file = file.path(output_data_dir, "CytoTRACE.rds"))

  phenot <- as.character(luad_epi$MP)
  names(phenot) <- rownames(luad_epi@meta.data)
  emb <- luad_epi@reductions[["umap.harmony"]]@cell.embeddings
  plotCytoTRACE(cyto_results, phenotype = phenot, emb = emb,
                outputDir = output_fig_dir)
  plotCytoGenes(cyto_results, numOfGenes = 30, outputDir = output_fig_dir)

  stem_score <- data.frame(
    cellnames = names(cyto_results[["CytoTRACE"]]),
    MP = luad_epi$MP,
    Stem_score = cyto_results[["CytoTRACE"]]
  )

  stem_score <- reorder_by_median(stem_score %>% rename(Score = Stem_score))
  stem_score <- stem_score %>% rename(Stem_score = Score)

  stem_colors <- c(
    "#3C5488FF", "#8491B4FF", "#4DBBD5FF", "#00A087FF",
    "#F39B7FFF", "#91D1C2FF", "#7E6148FF", "#DC0000FF", "#E64B35FF"
  )

  pdf(file.path(output_fig_dir, "Stem_Score.pdf"), width = 8, height = 4)
  print(plot_score_violin(stem_score %>% rename(Score = Stem_score), "Stem Score", stem_colors))
  dev.off()

  # MP proportion within clusters (pie charts)
  df <- data.frame(
    cellnames = colnames(cds),
    MP = cds@colData@listData[["MP"]],
    clusters = paste0("cluster_", cds@clusters@listData[["UMAP"]][["clusters"]])
  )
  df <- df[df$clusters != "cluster_7", ]
  df$MP <- factor(df$MP, levels = paste0("MP_", 1:9))
  df$clusters <- factor(df$clusters, levels = paste0("cluster_", 1:6))

  stage_colors <- c(
    "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
  )

  pdf(file.path(output_fig_dir, "MP_ratio_in_cluster.pdf"), width = 10, height = 8)
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))
  for (clust in paste0("cluster_", 1:6)) {
    mp_data <- df[df$clusters == clust, ]
    stage_counts <- table(mp_data$MP)
    pie(stage_counts, cex = 1,
        col = stage_colors,
        main = clust,
        labels = paste0(names(stage_counts), "\n", stage_counts, " (",
                        round(prop.table(stage_counts) * 100, 1), "%)"))
  }
  dev.off()

  message("CNV and stemness analysis complete.")
  return(list(cyto_trace = cyto_results))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
