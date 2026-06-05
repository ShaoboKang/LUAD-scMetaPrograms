# ============================================================================
# Script: 02_scissor_analysis
# Module: 08_clinical
# Description: Scissor analysis linking scRNA-seq MPs to bulk survival
#   phenotypes across TCGA+7 GEO cohorts.
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

#' Prepare and run Scissor for one dataset
#'
#' @param bulk_expr Bulk expression matrix.
#' @param sc_obj Seurat object.
#' @param phenotype Data frame with time and status columns.
#' @param dataset_name Name of the dataset.
#' @param output_data_dir Directory to save Scissor results.
#' @return Scissor result list.
run_scissor_dataset <- function(bulk_expr, sc_obj, phenotype, dataset_name, output_data_dir) {
  library(Scissor)

  save_file <- file.path(output_data_dir, paste0("Scissor_", dataset_name, "_survival.RData"))
  infos <- Scissor(bulk_expr, sc_obj, phenotype, alpha = 0.05,
                   family = "cox", Save_file = save_file)

  scissor_select <- rep("Background cells", ncol(sc_obj))
  names(scissor_select) <- colnames(sc_obj)
  scissor_select[infos$Scissor_pos] <- "Scissor + cells"
  scissor_select[infos$Scissor_neg] <- "Scissor - cells"

  message("Scissor for ", dataset_name, ": ")
  print(table(scissor_select))

  return(list(infos = infos, scissor_select = scissor_select))
}

# ---- Main Function ----

#' Run Scissor analysis across all bulk datasets
#'
#' @param sc_path Path to Seurat RDS.
#' @param bulk_path Path to TCGA+GEO bulk data.
#' @param scissor_func_path Path to Scissor utility script.
#' @param output_data_dir Directory for Scissor outputs.
#' @param output_fig_dir Directory for figures.
#' @return Seurat object with Scissor annotations.
run_analysis <- function(
    sc_path = file.path(RESULT_8_DIR, "results_data", "LUAD_Epi_assigned_MP_har.rds"),
    bulk_path = file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    scissor_func_path = file.path(CODE_DIR, "scripts", "08_clinical", "utils", "scissor_functions.R"),
    output_data_dir = file.path(RESULT_8_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(future)
  library(ggplot2)
  library(igraph)
  library(Seurat)
  library(Scissor)
  library(harmony)
  library(lemon)
  library(patchwork)

  if (exists("memory.limit", mode = "function")) {
    memory.limit(size = 900000000)
  }
  options(future.globals.maxSize = 50 * 1024^3)

  if (file.exists(scissor_func_path)) {
    source(scissor_func_path)
  } else {
    message("Scissor function script not found at: ", scissor_func_path)
  }

  message("Loading single-cell and bulk data...")
  luad_epi <- readRDS(sc_path)

  luad_epi$MP <- as.character(luad_epi$MP)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  luad_epi$MP <- unname(mp_map[luad_epi$MP])
  luad_epi$MP <- factor(luad_epi$MP, levels = paste0("MP_", 1:9))

  load(bulk_path)

  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  # Run Scissor for each dataset
  message("Running Scissor across datasets...")
  for (i in seq_along(expr_list)) {
    dataset_name <- names(expr_list)[i]
    message("Processing dataset: ", dataset_name)

    phenotype <- clin_list[[i]][, c("survival_time", "vital_status")]
    colnames(phenotype) <- c("time", "status")

    res <- run_scissor_dataset(expr_list[[i]], luad_epi, phenotype,
                               dataset_name, output_data_dir)

    luad_epi <- AddMetaData(luad_epi,
                            metadata = res$scissor_select,
                            col.name = paste0("Scissor_", dataset_name))
  }

  # Cell ratio bar plots
  ratio_plot <- list()
  col <- c("#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
           "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF")

  for (i in seq_along(expr_list)) {
    dataset_name <- names(expr_list)[i]
    meta_col <- paste0("Scissor_", dataset_name)

    cell_ratio <- as.data.frame(
      prop.table(table(luad_epi$MP, luad_epi@meta.data[[meta_col]]), margin = 2)
    )
    cell_ratio <- cell_ratio[-which(cell_ratio$Var2 == "Background cells"), ]

    p <- ggplot(cell_ratio) +
      geom_bar(aes(x = Var2, y = Freq, fill = Var1),
               stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
      theme_classic() +
      labs(title = dataset_name, x = "", y = "Ratio") +
      scale_fill_manual(values = col) +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

    ratio_plot[[dataset_name]] <- p
  }

  # UMAP and Scissor overlay plots
  p1 <- DimPlot(luad_epi, reduction = "umap.harmony",
                label = TRUE, repel = TRUE, label.size = 5, group.by = "MP") +
    scale_color_manual(values = col)

  p2 <- DimPlot(luad_epi, reduction = "umap.harmony", group.by = "Scissor_TCGA",
                cols = c("Scissor - cells" = "royalblue",
                         "Scissor + cells" = "indianred1",
                         "Background cells" = "grey"),
                pt.size = 1.2,
                order = c("Scissor - cells", "Scissor + cells"))

  pdf(file.path(output_fig_dir, "03_MP_scissor_TCGA.pdf"), width = 10, height = 4)
  print(p1 + p2)
  dev.off()

  pdf(file.path(output_fig_dir, "04_MP_scissor_ratio.pdf"), width = 22, height = 5)
  print(wrap_plots(ratio_plot, nrow = 1, ncol = 8) +
          plot_layout(guides = "collect") & theme(legend.position = "top"))
  dev.off()

  saveRDS(luad_epi, file = file.path(output_data_dir, "02_LUAD_Epi_assigned_MP_Scissor.rds"))

  message("Scissor analysis complete.")
  return(luad_epi)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
