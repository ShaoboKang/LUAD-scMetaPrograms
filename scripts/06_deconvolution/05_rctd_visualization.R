# ============================================================================
# Script: 05_rctd_visualization
# Module: 06_deconvolution
# Description:
#   Visualize RCTD spatial deconvolution results across Visium ST samples.
#   Generates spatial feature plots for selected MPs in cancer subregions.
#   Requires CARD_list.rds for subregion annotation (pre-computed).
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

#' Visualize RCTD spatial results
#'
#' @param rctd_path Path to RCTD results RDS.
#' @param card_path Path to CARD_list RDS (for subregion annotation).
#' @param output_fig_dir Directory for output figures.
#' @return Invisible NULL.
run_rctd_visualization <- function(
    rctd_path = file.path(RESULT_6_DIR, "results_data", "RCTD_visium_full.rds"),
    card_path = file.path(RESULT_6_DIR, "results_data", "CARD_list.rds"),
    output_fig_dir = file.path(RESULT_6_DIR, "results_figure")
) {
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(ggpubr)
  library(lemon)
  library(patchwork)
  library(RColorBrewer)

  if (!dir.exists(output_fig_dir)) {
    dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
  }

  message("Loading RCTD and CARD data...")
  CARD_list <- readRDS(card_path)
  RCTD <- readRDS(rctd_path)

  SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

  plot_sample_group <- function(samples, limits, output_file) {
    all_plots <- list()
    for (i in samples) {
      message("Plotting ", i, "...")
      normalize_weights <- as.data.frame(
        spacexr::normalize_weights(RCTD[[i]]@results$weights)
      )

      for (mp in c("MP_1", "MP_2", "MP_4", "MP_8")) {
        CARD_list[[i]][[mp]] <- normalize_weights[[mp]][
          match(colnames(CARD_list[[i]]), rownames(normalize_weights))
        ]
        CARD_list[[i]][[mp]][which(CARD_list[[i]]$subregion != "Cancer")] <- 0
        CARD_list[[i]][[mp]][is.na(CARD_list[[i]][[mp]])] <- 0
      }

      p1 <- SpatialPlot(CARD_list[[i]], group.by = "subregion",
                         label = TRUE, label.size = 3) +
        scale_color_npg() + scale_fill_npg() +
        ggplot2::theme(legend.position = "none")

      p2 <- SpatialFeaturePlot(CARD_list[[i]], features = "MP_1",
                                 ncol = 1, crop = TRUE) +
        scale_fill_gradientn(limits = c(0, limits[1]),
                             colours = SpatialColors(n = 100))
      p3 <- SpatialFeaturePlot(CARD_list[[i]], features = "MP_2",
                                 ncol = 1, crop = TRUE) +
        scale_fill_gradientn(limits = c(0, limits[2]),
                             colours = SpatialColors(n = 100))
      p4 <- SpatialFeaturePlot(CARD_list[[i]], features = "MP_4",
                                 ncol = 1, crop = TRUE) +
        scale_fill_gradientn(limits = c(0, limits[3]),
                             colours = SpatialColors(n = 100))
      p5 <- SpatialFeaturePlot(CARD_list[[i]], features = "MP_8",
                                 ncol = 1, crop = TRUE) +
        scale_fill_gradientn(limits = c(0, limits[4]),
                             colours = SpatialColors(n = 100))

      sample_plots <- (p1 + p2 + p3 + p4 + p5) +
        plot_layout(widths = c(1, 1, 1, 1, 1))
      all_plots[[i]] <- sample_plots
    }

    final_plot <- wrap_plots(all_plots, nrow = length(all_plots))
    pdf(output_file, height = 20, width = 12)
    print(final_plot)
    dev.off()
    message("Saved: ", output_file)
  }

  # Group 1
  plot_sample_group(
    c("seurat_ST8_AIS", "seurat_ST3_MIA", "seurat_ST1_IAC"),
    limits = c(0.5, 0.6, 1, 0.6),
    file.path(output_fig_dir, "all_samples_combined_ratio.pdf")
  )

  # Group 2
  plot_sample_group(
    c("seurat_ST5_AIS", "seurat_ST6_MIA", "seurat_ST2_IAC"),
    limits = c(0.2, 0.4, 1, 0.6),
    file.path(output_fig_dir, "all_samples_combined_ratio2.pdf")
  )

  message("[INFO] RCTD visualization complete.")
  invisible(NULL)
}

# ---- Entry Point ----
if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_rctd_visualization()
} else {
  message("[INFO] Script loaded. Call run_rctd_visualization() to execute.")
}
