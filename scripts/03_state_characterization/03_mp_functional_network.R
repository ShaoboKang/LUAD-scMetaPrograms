# ============================================================================
# Script: 03_mp_functional_network
# Module: 03_state_characterization
# Description: Functional enrichment (GO-BP) comparison across MP gene sets
#   using compareCluster and emapplot visualization.
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

# ---- Main Function ----

#' Run functional enrichment comparison across MPs
#'
#' @param mp_list_path Path to MP DEG list RDS.
#' @param tf_path Path to top5 TF RDS.
#' @param output_fig_dir Directory for figure output.
#' @return compareCluster result object.
run_analysis <- function(
    mp_list_path = file.path(RESULT_3_DIR, "results_data", "01_tumor_DEG_list.rds"),
    tf_path = file.path(RESULT_4_DIR, "results_data", "01_top5_TF.rds"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
  library(dplyr)
  library(GOSemSim)
  library(ggforce)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading MP gene list and TF names...")
  mp_list <- readRDS(file = mp_list_path)

  tf_names <- readRDS(tf_path)
  tf_names <- unique(gsub("\\(.*\\)", "", tf_names$Topic))

  message("Running compareCluster GO enrichment...")
  result <- compareCluster(
    geneCluster = mp_list,
    fun = "enrichGO",
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 0.05,
    maxGSSize = 2000,
    ont = "BP"
  )

  result@compareClusterResult <- result@compareClusterResult %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()

  d <- godata("org.Hs.eg.db", ont = "BP")
  result2 <- pairwise_termsim(result, method = "JC", semData = d)

  p_final <- emapplot(
    result2,
    showCategory = 200,
    pie = "equal",
    layout = "kk",
    size_category = 10,
    size_edge = 0.5,
    node_label = "category",
    min_edge = 0.5,
    color_edge = "black",
    group = TRUE,
    group_style = "polygon",
    clusterFunction = cluster::clara,
    nCluster = 8
  )

  a <- p_final$data
  message("Unique cluster labels: ", length(unique(a$label)))

  my_colors <- c(
    "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
  )
  names(my_colors) <- paste0("MP_", 1:9)

  label_positions <- a %>%
    group_by(color2) %>%
    summarise(
      x = mean(x),
      y = mean(y),
      .groups = "drop"
    ) %>%
    mutate(cluster_name = color2)
  message("Unique label positions: ", length(unique(label_positions$cluster_name)))

  # Re-run emapplot without node labels for custom annotation
  p_final <- emapplot(
    result2,
    showCategory = 200,
    pie = "equal",
    layout = "kk",
    size_category = 10,
    size_edge = 0.2,
    node_label = "category",
    min_edge = 0.5,
    color_edge = "grey",
    group = FALSE
  )

  aa <- result2@compareClusterResult$Count[
    match(p_final$layers[[2]]$data$label, result2@compareClusterResult$Description)
  ]
  aa1 <- c(max(log(aa)), min(log(aa)))

  point_size <- 0.05
  p_final@layers[[2]]$data$pathway_size <- log(aa) / (aa1[1] - aa1[2]) * point_size
  message("Unique pathways: ", length(unique(p_final$layers[[2]]$data$name)))

  col <- rep("#8491B4FF", 35)

  plot <- p_final +
    ggforce::geom_mark_hull(
      data = a,
      aes(x, y, group = color2, colour = color2),
      alpha = 0.1, linewidth = 1, linetype = "dashed",
      expand = unit(2.5, "pt"),
      concavity = 4
    ) +
    scale_color_manual(values = col) +
    geom_text(
      data = label_positions,
      aes(x = x, y = y, label = cluster_name),
      fontface = "bold",
      size = 8
    ) +
    scale_fill_manual(values = my_colors) +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "#F8FCFE")
    )

  pdf(file.path(output_fig_dir, "MP_function_legend1.pdf"),
      width = 15, height = 15, onefile = FALSE)
  print(plot)
  dev.off()

  message("Functional differences analysis complete.")
  return(result2)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
