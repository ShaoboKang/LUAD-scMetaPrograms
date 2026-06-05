#' VPS72 Spatial Transcriptomics Analysis
# Script: 03_vps72_spatial
#'
#' @description
#' Visualizes VPS72 expression across spatial transcriptomics subregions
#' using CARD-deconvoluted Seurat objects.
#'
#' @return Invisibly returns NULL. Writes figures to RESULT_4_DIR.
#'
#' @examples
#' \dontrun{
#' run_vps72_spatial()
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
  library(ggsci)
  library(ggpubr)
  library(lemon)
  library(patchwork)
})

#' Plot VPS72 Spatial Expression for One Sample
#'
#' @param card_obj Seurat object from CARD deconvolution.
#' @param sample_name Sample identifier string.
#' @return Patchwork plot of spatial subregion, feature plot, and boxplot.
plot_vps72_per_sample <- function(card_obj, sample_name) {
  p0 <- SpatialPlot(card_obj, group.by = "subregion", label = TRUE, label.size = 5) +
    scale_color_npg() +
    scale_fill_npg() +
    ggplot2::theme(legend.position = "none")

  p1 <- SpatialFeaturePlot(card_obj, features = c("VPS72"))

  combined_df <- data.frame(
    VPS72 = card_obj@assays[["SCT"]]@data["VPS72", ],
    subregion = card_obj$subregion
  )

  combined_df$subregion <- factor(combined_df$subregion, levels = unique(combined_df$subregion))
  subregion <- levels(factor(combined_df$subregion))
  comp <- combn(subregion, 2)
  my_comparisons <- list()
  for (j in seq_len(ncol(comp))) {
    my_comparisons[[j]] <- comp[, j]
  }

  p2 <- ggboxplot(
    combined_df, "subregion", "VPS72",
    fill = "subregion", outlier.shape = NA, title = sample_name
  ) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.title.x = element_blank()) +
    stat_compare_means(
      method = "t.test", hide.ns = FALSE,
      comparisons = my_comparisons, label = "p.adj.format",
      p.adjust.method = "BH", vjust = 0.02, bracket.size = 0.6
    ) +
    scale_fill_manual(values = c("#E64B35FF", "#79B99D", "#4DBBD5FF", "#fce38a"))

  p0 + p1 + p2
}

#' Run VPS72 Spatial Analysis
#'
#' @return Invisibly returns NULL.
run_vps72_spatial <- function() {
  message("[03_vps72_spatial] Starting spatial analysis...")

  card_list <- readRDS(file = file.path(RESULT_6_DIR, "results_data", "CARD_list.rds"))

  samples <- c("seurat_ST8_AIS", "seurat_ST3_MIA", "seurat_ST1_IAC")

  for (i in samples) {
    p_combined <- plot_vps72_per_sample(card_list[[i]], i)
    print(p_combined)
  }

  pdf(file = file.path(RESULT_4_DIR, "results_figure", "09_VPS72_exp_in_ST.pdf"), width = 12, height = 10)
  print(p_combined)
  dev.off()

  message("[03_vps72_spatial] Completed successfully.")
  invisible(NULL)
}

if (!interactive()) {
  run_vps72_spatial()
}
