#' Spatial Expression of Ligand-Receptor Pairs
# Script: 08_lr_spatial_expression
#'
#' Analyzes spatial transcriptomics data to visualize MDK grouping and
#' gene-pair co-expression patterns. Generates spatial plots and Venn diagrams
#' with hypergeometric test p-values for each sample.
#'
#' @return None. PDF figures are written to the result figure directory.
#' @author Script refactored from original analysis pipeline.
#' @export
run_lr_spatial_expression <- function() {
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
  .libPaths(c("~/R/SeuratV4", .libPaths()))
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(ggpubr)
  library(lemon)
  library(patchwork)
  library(ggvenn)

  CARD_list <- readRDS(file = file.path(RESULT_6_DIR, "results_data", "CARD_list.rds"))

  # Define gene pairs to analyze
  gene_pairs <- list(
    c("ANGPTL4", "CDH5"),
    c("ANGPTL4", "CDH11"),
    c("SCGB3A2", "MARCO"),
    c("PLAU", "PLAUR")
  )

  for (i in names(CARD_list)) {
    # Get all required genes
    all_genes <- unique(c("MDK", unlist(gene_pairs)))

    # Fetch gene expression data
    mecom_expr <- FetchData(CARD_list[[i]], vars = all_genes, layer = "data")

    # Create MDK grouping
    mecom_expr <- mecom_expr %>%
      mutate(
        MDK_group = ifelse(MDK > 0, "MDK+", "MDK-")
      )

    # Add MDK group to metadata
    CARD_list[[i]]$MDK_group <- factor(mecom_expr$MDK_group, levels = c("MDK+", "MDK-"))

    # Plot MDK grouping
    p_MDK <- SpatialPlot(CARD_list[[i]], group.by = "MDK_group") +
      scale_fill_manual(values = c("#DC0000FF", "grey")) +
      theme(legend.position = "right")

    # Create grouping and plots for each gene pair
    plots_list <- list()
    plots_list[[1]] <- p_MDK

    plots_list[[2]] <- SpatialPlot(CARD_list[[i]], group.by = "subregion") +
      scale_color_npg() +
      scale_fill_npg() +
      theme(legend.position = "right")

    plot_index <- 3

    for (j in seq_along(gene_pairs)) {
      gene_pair <- gene_pairs[[j]]
      gene1 <- gene_pair[1]
      gene2 <- gene_pair[2]

      # Create grouping based on gene pair
      group_col_name <- paste0(gene1, "_", gene2, "_group")

      mecom_expr <- mecom_expr %>%
        mutate(
          !!group_col_name := case_when(
            MDK_group == "MDK-" ~ "MDK-",
            !!sym(gene1) > 0 & !!sym(gene2) > 0 ~ paste0(gene1, "+", gene2, "+"),
            !!sym(gene1) > 0 & !!sym(gene2) == 0 ~ paste0(gene1, "+"),
            !!sym(gene1) == 0 & !!sym(gene2) > 0 ~ paste0(gene2, "+"),
            !!sym(gene1) == 0 & !!sym(gene2) == 0 ~ paste0(gene1, "-", gene2, "-"),
            TRUE ~ "MDK-"
          )
        )

      # Convert grouping to factor with defined order
      group_levels <- c(
        paste0(gene1, "+", gene2, "+"),
        paste0(gene1, "+"),
        paste0(gene2, "+"),
        paste0(gene1, "-", gene2, "-"),
        "MDK-"
      )

      mecom_expr[[group_col_name]] <- factor(mecom_expr[[group_col_name]], levels = group_levels)

      # Add to Seurat object
      CARD_list[[i]][[group_col_name]] <- mecom_expr[[group_col_name]]

      # Color scheme for this gene pair
      group_colors <- c("#DC0000FF", "#3C5488FF", "#8491B4FF", "#91D1C2FF", "grey")

      # Plot grouping spatial map
      p_group <- SpatialPlot(CARD_list[[i]], group.by = group_col_name) +
        scale_fill_manual(values = group_colors) +
        theme(legend.position = "right")

      # Create Venn diagram (MDK+ cells only)
      mdk_pos_cells <- rownames(mecom_expr)[mecom_expr$MDK_group == "MDK+"]
      mecom_expr_mdk_pos <- mecom_expr[mdk_pos_cells, ]

      venn_list <- setNames(
        list(
          rownames(mecom_expr_mdk_pos)[mecom_expr_mdk_pos[[group_col_name]] %in%
                                         c(paste0(gene1, "+"), paste0(gene1, "+", gene2, "+"))],
          rownames(mecom_expr_mdk_pos)[mecom_expr_mdk_pos[[group_col_name]] %in%
                                         c(paste0(gene2, "+"), paste0(gene1, "+", gene2, "+"))]
        ),
        c(paste0(gene1, "+"), paste0(gene2, "+"))
      )

      # Hypergeometric test for overlap significance
      N <- nrow(mecom_expr_mdk_pos)
      K <- length(venn_list[[1]])
      n <- length(venn_list[[2]])
      k <- length(intersect(venn_list[[1]], venn_list[[2]]))
      p_value <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
      p_value <- sprintf("%.3f", p_value)

      # Plot Venn diagram
      p_venn <- ggvenn(
        venn_list,
        fill_color = c("#3C5488FF", "#8491B4FF"),
        stroke_linetype = "longdash"
      ) +
        labs(caption = paste("p =", p_value)) +
        theme(
          plot.caption = element_text(
            hjust = 0.5,
            size = 10,
            face = "italic",
            margin = margin(t = 10)
          )
        )

      # Store plots
      plots_list[[plot_index]] <- p_group
      plots_list[[plot_index + 1]] <- p_venn
      plot_index <- plot_index + 2
    }

    # Final layout: 5 columns x 2 rows
    final_layout <- wrap_plots(
      list(
        plots_list[[1]], plots_list[[3]], plots_list[[5]], plots_list[[7]], plots_list[[9]],
        plots_list[[2]], plots_list[[4]], plots_list[[6]], plots_list[[8]], plots_list[[10]]
      ),
      ncol = 5,
      nrow = 2,
      widths = rep(1, 5),
      heights = c(1, 2)
    )

    # Save as PDF
    output_dir <- file.path(RESULT_5_DIR, "results_figure")
    pdf_file <- file.path(output_dir, paste0("Sample_", i, "_MDK_gene_pairs_analysis.pdf"))

    pdf(pdf_file, width = 22, height = 7)
    print(final_layout)
    dev.off()

    cat("Sample", i, "analysis completed. Plots saved to:", pdf_file, "\n")
  }
}

if (!interactive()) {
  run_lr_spatial_expression()
}
