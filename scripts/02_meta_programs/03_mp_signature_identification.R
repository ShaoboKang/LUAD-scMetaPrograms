# ============================================================================
# Script: 03_mp_signature_identification
# Module: 02_meta_programs
# Description:
#   Load NMF results, select robust NMF programs, cluster into Meta-Programs (MPs),
#   compute MP abundances, assign MP labels to tumor cells, and generate visualizations.
#
# Prerequisites:
#   NMF results must be generated first by running:
#     02_nmf_per_sample.R (server-side NMF computation)
#   This produces:
#     Result_2/results_data/NMF_results/<sample>_rank4_9_nruns10.RDS
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
library(parallel)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

# ---- Helper Functions ----

#' Load NMF result files and extract W basis matrices.
#' @param nmf_dir Directory containing NMF result RDS files.
#' @return Named list of W basis matrices.
load_nmf_results <- function(nmf_dir) {
  f_sam <- dir(nmf_dir, pattern = "_rank4_9_nruns10\\.RDS")
  genes_nmf_w_basis_rank <- lapply(f_sam, function(sam) {
    nmf_res <- readRDS(file.path(nmf_dir, sam))
    valid_ranks <- as.numeric(names(which(is.na(nmf_res[["fit"]]) == FALSE)))
    w_basis <- sapply(valid_ranks, function(n_rank) {
      w <- NMF::basis(nmf_res$fit[[as.character(n_rank)]])
      colnames(w) <- paste(sam, n_rank, 1:n_rank, sep = ".")
      return(w)
    })
    w_basis <- do.call(cbind, w_basis)
    w_basis
  })
  names(genes_nmf_w_basis_rank) <- f_sam
  genes_nmf_w_basis_rank
}

#' Extract top genes for each NMF program.
#' @param genes_nmf_w_basis_rank List of W basis matrices.
#' @param program_size Number of top genes per program. Default 50.
#' @return List of top gene matrices.
extract_nmf_programs <- function(genes_nmf_w_basis_rank, program_size = 50) {
  nmf_programs <- lapply(genes_nmf_w_basis_rank, function(x) {
    apply(x, 2, function(y) names(sort(y, decreasing = TRUE))[1:program_size])
  })
  nmf_programs <- lapply(nmf_programs, toupper)
  nmf_programs
}

#' Run robust NMF program filtering.
#' @param nmf_programs List of NMF program matrices.
#' @param intra_min Minimum intra-sample overlap. Default 35.
#' @param intra_max Maximum intra-sample overlap for redundancy removal. Default 10.
#' @param inter_min Minimum inter-sample overlap. Default 10.
#' @return Character vector of selected program names.
filter_robust_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_min = 10) {
  robust_nmf_programs(nmf_programs, intra_min = intra_min, intra_max = intra_max,
                      inter_filter = TRUE, inter_min = inter_min)
}

#' Build MP clusters from NMF program similarity.
#' @param nmf_programs Matrix of NMF programs (genes x programs).
#' @param program_size Program size for border gene handling. Default 50.
#' @param min_intersect_initial Initial intersection threshold. Default 10.
#' @param min_intersect_cluster Cluster addition threshold. Default 10.
#' @param min_group_size Minimum group size. Default 5.
#' @param genes_nmf_w_basis_rank Original W basis list for score-based tie-breaking.
#' @return List with Cluster_list and MP_list.
build_mp_clusters <- function(nmf_programs, program_size = 50,
                              min_intersect_initial = 10, min_intersect_cluster = 10,
                              min_group_size = 5, genes_nmf_w_basis_rank) {
  nmf_intersect <- parallel_intersection(nmf_programs, program_size)
  nmf_intersect_hc <- hclust(as.dist(program_size - nmf_intersect), method = "average")
  nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
  nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

  sorted_intersection <- sort(apply(nmf_intersect, 2, function(x) (length(which(x >= min_intersect_initial)) - 1)), decreasing = TRUE)

  cluster_list <- list()
  mp_list <- list()
  k <- 1
  curr_cluster <- c()
  nmf_intersect_original <- nmf_intersect

  while (sorted_intersection[1] > min_group_size) {
    curr_cluster <- c(curr_cluster, names(sorted_intersection[1]))
    genes_mp <- nmf_programs[, names(sorted_intersection[1])]
    nmf_programs <- nmf_programs[, -match(names(sorted_intersection[1]), colnames(nmf_programs))]
    intersection_with_genes_mp <- sort(apply(nmf_programs, 2, function(x) length(intersect(genes_mp, x))), decreasing = TRUE)
    nmf_history <- genes_mp

    while (intersection_with_genes_mp[1] >= min_intersect_cluster) {
      curr_cluster <- c(curr_cluster, names(intersection_with_genes_mp)[1])
      genes_mp_temp <- sort(table(c(nmf_history, nmf_programs[, names(intersection_with_genes_mp)[1]])), decreasing = TRUE)
      genes_at_border <- genes_mp_temp[which(genes_mp_temp == genes_mp_temp[program_size])]

      if (length(genes_at_border) > 1) {
        genes_curr_nmf_score <- c()
        for (i in curr_cluster) {
          curr_study <- paste(strsplit(i, "[.]")[[1]][1:which(strsplit(i, "[.]")[[1]] == "RDS")], collapse = ".")
          q <- genes_nmf_w_basis_rank[[curr_study]][
            match(names(genes_at_border), toupper(rownames(genes_nmf_w_basis_rank[[curr_study]])))[!is.na(match(names(genes_at_border), toupper(rownames(genes_nmf_w_basis_rank[[curr_study]]))))],
            i
          ]
          names(q) <- names(genes_at_border[!is.na(match(names(genes_at_border), toupper(rownames(genes_nmf_w_basis_rank[[curr_study]]))))])
          genes_curr_nmf_score <- c(genes_curr_nmf_score, q)
        }
        genes_curr_nmf_score_sort <- sort(genes_curr_nmf_score, decreasing = TRUE)
        genes_curr_nmf_score_sort <- genes_curr_nmf_score_sort[unique(names(genes_curr_nmf_score_sort))]
        genes_mp_temp <- c(names(genes_mp_temp[which(genes_mp_temp > genes_mp_temp[program_size])]), names(genes_curr_nmf_score_sort))
      } else {
        genes_mp_temp <- names(genes_mp_temp)[1:program_size]
      }

      nmf_history <- c(nmf_history, nmf_programs[, names(intersection_with_genes_mp)[1]])
      genes_mp <- genes_mp_temp[1:program_size]
      nmf_programs <- nmf_programs[, -match(names(intersection_with_genes_mp)[1], colnames(nmf_programs))]
      intersection_with_genes_mp <- sort(apply(nmf_programs, 2, function(x) length(intersect(genes_mp, x))), decreasing = TRUE)
    }

    cluster_list[[paste0("Cluster_", k)]] <- curr_cluster
    mp_list[[paste0("MP_", k)]] <- genes_mp
    k <- k + 1

    nmf_intersect <- nmf_intersect[-match(curr_cluster, rownames(nmf_intersect)), -match(curr_cluster, colnames(nmf_intersect))]
    sorted_intersection <- sort(apply(nmf_intersect, 2, function(x) (length(which(x >= min_intersect_initial)) - 1)), decreasing = TRUE)
    curr_cluster <- c()
    message("Remaining programs: ", dim(nmf_intersect)[2])
  }

  list(Cluster_list = cluster_list, MP_list = mp_list, nmf_intersect_original = nmf_intersect_original)
}

#' Compute pairwise intersection matrix in parallel.
#' @param nmf_programs Matrix of programs.
#' @param program_size Program size for intersection calculation. Default 50.
#' @return Intersection matrix.
parallel_intersection <- function(nmf_programs, program_size = 50) {
  cl <- makeCluster(8, type = "FORK")
  nmf_intersect <- parApply(cl = cl, nmf_programs, 2, function(x) {
    apply(nmf_programs, 2, function(y) length(intersect(x, y)))
  })
  stopCluster(cl)
  nmf_intersect
}

#' Plot Jaccard similarity heatmap of NMF programs.
#' @param nmf_intersect_original Original intersection matrix.
#' @param cluster_list Cluster list.
#' @param output_path Output PDF path.
#' @param custom_magma Custom color palette.
plot_jaccard_heatmap <- function(nmf_intersect_original, cluster_list, output_path, custom_magma) {
  new_order <- c("Cluster_1", "Cluster_8", "Cluster_2", "Cluster_7", "Cluster_4",
                 "Cluster_6", "Cluster_5", "Cluster_3", "Cluster_9")
  cluster_list <- cluster_list[new_order]

  inds_sorted <- c()
  for (j in 1:length(cluster_list)) {
    inds_sorted <- c(inds_sorted, match(cluster_list[[j]], colnames(nmf_intersect_original)))
  }
  inds_new <- c(inds_sorted, which(is.na(match(1:dim(nmf_intersect_original)[2], inds_sorted))))

  nmf_intersect_meltI_new <- reshape2::melt(nmf_intersect_original[inds_sorted, rev(inds_sorted)])

  p <- ggplot(data = nmf_intersect_meltI_new, aes(x = Var1, y = Var2, fill = 100 * value / (100 - value), color = 100 * value / (100 - value))) +
    geom_tile() +
    scale_color_gradient2(limits = c(2, 25), low = custom_magma[1:111], mid = custom_magma[112:222],
                          high = custom_magma[223:333], midpoint = 13.5, oob = squish, name = "Similarity\n(Jaccard index)") +
    scale_fill_gradient2(limits = c(2, 25), low = custom_magma[1:111], mid = custom_magma[112:222],
                         high = custom_magma[223:333], midpoint = 13.5, oob = squish, name = "Similarity\n(Jaccard index)") +
    theme(axis.ticks = element_blank(), panel.border = element_rect(fill = FALSE),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text = element_text(size = 11), axis.title = element_text(size = 12),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10),
          legend.text.align = 0.5, legend.justification = "bottom") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

  pdf(file = output_path, width = 8, height = 6)
  print(p)
  dev.off()
  message("Saved Jaccard heatmap to: ", output_path)
}

#' Calculate MP abundance across stages.
#' @param cluster_list Cluster list.
#' @param mp_list MP list.
#' @param nmf_programs_original Original NMF programs matrix.
#' @param luad_epi Seurat object for sample metadata.
#' @return Data frame of abundance statistics.
calculate_mp_abundance <- function(cluster_list, mp_list, nmf_programs_original, luad_epi) {
  programs_info <- data.frame(names = colnames(nmf_programs_original))
  programs_info$Sample <- sub("__.+", "", programs_info$names)
  programs_info$stage <- luad_epi$stage[match(programs_info$Sample, luad_epi$Sample)]
  cancer_list <- split(programs_info$names, programs_info$stage)

  abundance <- plyr::ldply(1:length(cluster_list), function(i) {
    mmp <- names(mp_list)[i]
    plyr::ldply(names(cancer_list), function(cancer) {
      x <- cluster_list[[i]]
      y <- cancer_list[[cancer]]
      observed <- length(intersect(x, y))
      expected <- length(x) * length(y) / length(unlist(cancer_list))
      a <- log2((observed + 1) / (expected + 1))
      p <- phyper(observed - 1, length(y), length(unlist(cancer_list)) - length(y), length(x), lower.tail = FALSE)
      data.frame(Cluster = mmp, Cancer = cancer, observed = observed, expected = expected,
                 MP_size = length(x), cancer_size = length(y), A = a, p_value = p)
    })
  })

  abundance$p_adj <- p.adjust(abundance$p_value, method = "bonferroni")
  abundance <- mutate(abundance,
    classification = ifelse((observed > 10) | (A > 1),
      ifelse(p_adj < 0.05, "High_significant", "High"),
      ifelse(((observed <= 10) & (observed >= 2)) | ((A <= 1) & (A > 0)), "Medium",
        ifelse((observed == 1) & ((A <= 0) & (A > -1.5)), "Low",
          ifelse(observed == 0, "Absent", "Absent")))))
  abundance$classification <- factor(abundance$classification,
                                     levels = c("Absent", "Low", "Medium", "High", "High_significant"))
  abundance$Cluster <- factor(abundance$Cluster, levels = paste0("MP_", 1:9))
  abundance$Cancer <- factor(abundance$Cancer, levels = c("Early", "Advance"))
  abundance
}

#' Plot MP abundance heatmap.
#' @param abundance Abundance data frame.
#' @param output_path Output PDF path.
plot_mp_abundance <- function(abundance, output_path) {
  my_cols <- c("white", viridis::magma(4, direction = -1, begin = 0.18, end = 1))
  names(my_cols) <- c("Absent", "Low", "Medium", "High", "High_significant")

  plot_data <- reshape2::dcast(abundance, Cluster ~ Cancer, value.var = "classification")
  rownames(plot_data) <- plot_data$Cluster
  tmp_text <- grid::gpar(fontfamily = "sans", fontsize = 8)

  pdf(file = output_path, width = 8, height = 6)
  print(Heatmap(plot_data[rev(1:nrow(plot_data)), -1], col = my_cols, border = TRUE,
                row_names_gp = tmp_text, column_names_gp = tmp_text, row_names_side = "left",
                row_title_gp = tmp_text, column_title_gp = tmp_text,
                height = nrow(plot_data) * unit(0.3, "cm"),
                width = ncol(plot_data) * unit(0.3, "cm"),
                heatmap_legend_param = list(border = "black"),
                name = " ", cluster_rows = FALSE, cluster_columns = FALSE))
  dev.off()
  message("Saved MP abundance heatmap to: ", output_path)
}

#' Plot MP stage distribution barplot.
#' @param cluster_list Cluster list.
#' @param luad_epi Seurat object.
#' @param output_path Output PDF path.
plot_mp_stage_distribution <- function(cluster_list, luad_epi, output_path) {
  cluster_list_abund <- lapply(cluster_list, function(x) {
    data.frame(row.names = x, sample = sub("__.+", "", x),
               stage = luad_epi$stage[match(sub("__.+", "", x), luad_epi$Sample)])
  })
  for (i in 1:length(cluster_list_abund)) {
    cluster_list_abund[[i]]$MP <- paste0("MP", "_", i)
  }
  cluster_list_abund <- do.call(rbind, cluster_list_abund)
  cluster_list_abund$stage <- factor(cluster_list_abund$stage, levels = c("Early", "Advance"))
  count_table <- table(cluster_list_abund$stage, cluster_list_abund$MP) %>% as.data.frame()
  names(count_table) <- c("stage", "MP", "count")

  p <- ggplot(count_table, aes(x = MP, y = count, fill = stage)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = count), position = position_dodge(0.8), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("#1F77B4", "#D62728")) +
    labs(x = "MP", y = "Counts") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  pdf(file = output_path, width = 10, height = 6)
  print(p)
  dev.off()
  message("Saved MP stage distribution to: ", output_path)
}

#' Assign MP labels to cells based on signature scores.
#' @param expr_list List of expression matrices per sample.
#' @param mp_list MP gene lists.
#' @param min_genes Minimum conserved genes fraction. Default 0.7.
#' @param min_score Minimum score threshold for assignment. Default 0.8.
#' @return List of MP label vectors per sample.
assign_mp_labels <- function(expr_list, mp_list, min_genes = 0.7, min_score = 0.8) {
  mp_scores_per_sample <- lapply(expr_list, function(x) {
    sigScores(x, sigs = mp_list, conserved.genes = min_genes)
  })

  assign_mp_per_cell <- lapply(mp_scores_per_sample, function(x) {
    apply(x, 1, function(y) {
      max_val <- max(y, na.rm = TRUE)
      max_idx <- which.max(y)
      if (max_val < min_score) {
        return("un_assigned")
      } else {
        return(colnames(x)[max_idx])
      }
    })
  })

  list(mp_scores = mp_scores_per_sample, labels = assign_mp_per_cell)
}

#' Plot stacked barplot of MP proportions per sample.
#' @param df Data frame with MP and sample columns.
#' @param luad_epi Seurat object for stage metadata.
#' @param output_path Output PDF path.
#' @param col Color palette.
#' @param include_unassigned Logical; include unassigned cells.
plot_mp_sample_barplot <- function(df, luad_epi, output_path, col, include_unassigned = TRUE) {
  count_table <- table(df$MP, df$sample) %>% as.data.frame()
  names(count_table) <- c("MP", "sample", "count")
  count_table <- count_table %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

  df_wide <- count_table %>%
    pivot_wider(id_cols = sample, names_from = MP, values_from = proportion)
  df_wide$stage <- luad_epi$stage[match(df_wide$sample, luad_epi$Sample)]

  if (include_unassigned) {
    df_wide <- df_wide %>%
      arrange(stage, un_assigned, MP_1, MP_2, MP_3, MP_4, MP_5, MP_6, MP_7, MP_8, MP_9)
  } else {
    df_wide <- df_wide %>%
      arrange(stage, MP_1, MP_2, MP_3, MP_4, MP_5, MP_6, MP_7, MP_8, MP_9)
  }

  sample_order <- df_wide %>% pull(sample)
  count_table <- count_table %>%
    mutate(sample = factor(sample, levels = sample_order))

  p <- ggplot(count_table, aes(x = sample, y = proportion, fill = factor(MP))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = col, name = "MP") +
    labs(x = "Sample ID", y = "Proportion", title = "Proportion of MP Tumor Cells by Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position = "top")

  pdf(file = output_path, width = 10, height = 6)
  print(p)
  dev.off()
  message("Saved MP sample barplot to: ", output_path)
}

#' Plot pie charts of MP cells by stage.
#' @param df Data frame with stage and MP_anno columns.
#' @param output_path Output PDF path.
plot_mp_stage_pies <- function(df, output_path) {
  stage_colors <- c("Early" = "#1F77B4", "Advance" = "#D62728")
  pdf(file = output_path, width = 8, height = 8)
  par(mfrow = c(3, 3), mar = c(2, 2, 3, 2))
  for (mp in paste0("MP_", 1:9)) {
    mp_data <- df[df$MP_anno == mp, ]
    stage_counts <- table(mp_data$stage)
    pie(stage_counts, cex = 1, col = stage_colors[names(stage_counts)],
        main = mp,
        labels = paste0(names(stage_counts), "\n", stage_counts, " (",
                        round(prop.table(stage_counts) * 100, 1), "%)"))
  }
  dev.off()
  message("Saved MP stage pie charts to: ", output_path)
}

# ---- Main Function ----

#' Run MP signature identification pipeline.
#' @param nmf_dir Directory with NMF result files.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return List of results including Cluster_list, MP_list, and assigned labels.
run_analysis <- function(
    nmf_dir = file.path(RESULT_2_DIR, "results_data", "NMF_results"),
    output_data_dir = file.path(RESULT_2_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  # Load utility functions
  source(file.path(CODE_DIR, "scripts", "02_meta_programs", "utils", "00_tools_robust_nmf_programs.R"))

  # Load NMF results
  message("Loading NMF results...")
  genes_nmf_w_basis_rank <- load_nmf_results(nmf_dir)

  # Extract top programs
  program_size <- 50
  intra_min_parameter <- 35
  intra_max_parameter <- 10
  inter_min_parameter <- 10

  nmf_programs <- extract_nmf_programs(genes_nmf_w_basis_rank, program_size = program_size)

  # Filter robust programs
  message("Filtering robust NMF programs...")
  nmf_filter_ccle <- filter_robust_programs(nmf_programs, intra_min = intra_min_parameter,
                                            intra_max = intra_max_parameter, inter_min = inter_min_parameter)
  nmf_programs <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle), drop = FALSE])
  nmf_programs <- do.call(cbind, nmf_programs)
  nmf_programs_original <- nmf_programs

  # Build MP clusters
  message("Clustering NMF programs into MPs...")
  cluster_res <- build_mp_clusters(nmf_programs, program_size = program_size,
                                   genes_nmf_w_basis_rank = genes_nmf_w_basis_rank)
  cluster_list <- cluster_res$Cluster_list
  mp_list <- cluster_res$MP_list
  nmf_intersect_original <- cluster_res$nmf_intersect_original

  # Reorder clusters
  new_order <- c("Cluster_1", "Cluster_8", "Cluster_2", "Cluster_7", "Cluster_4",
                 "Cluster_6", "Cluster_5", "Cluster_3", "Cluster_9")
  cluster_list <- cluster_list[new_order]
  mp_list <- mp_list[new_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))

  # Plot Jaccard heatmap
  plot_jaccard_heatmap(nmf_intersect_original, cluster_list,
                       file.path(output_fig_dir, "02_MP_cluster.pdf"),
                       custom_magma)

  # Save clusters
  save(Cluster_list = cluster_list, MP_list = mp_list, nmf_programs_original,
       nmf_intersect_original, file = file.path(output_data_dir, "02_Malig_final_MP_top50.RData"))

  # Load LUAD data for abundance and assignment
  luad_epi <- readRDS(file.path(output_data_dir, "01_LUAD_epi_anno_final_NMF_input.rds"))

  # MP abundance
  message("Calculating MP abundances...")
  abundance <- calculate_mp_abundance(cluster_list, mp_list, nmf_programs_original, luad_epi)
  plot_mp_abundance(abundance, file.path(output_fig_dir, "03_MP_cluster_abundance.pdf"))
  plot_mp_stage_distribution(cluster_list, luad_epi, file.path(output_fig_dir, "04_MP_cluster_stage.pdf"))

  # Load scoring utilities
  source(file.path(CODE_DIR, "scripts", "02_meta_programs", "utils", "matrix-ops.R"))
  source(file.path(CODE_DIR, "scripts", "02_meta_programs", "utils", "sigScores.R"))
  source(file.path(CODE_DIR, "scripts", "02_meta_programs", "utils", "bin.R"))
  source(file.path(CODE_DIR, "scripts", "02_meta_programs", "utils", "utils-bin.R"))

  # Assign MP labels per sample
  message("Assigning MP labels to cells...")
  min_genes <- 35
  min_score <- 0.8

  my_study <- SplitObject(luad_epi, split.by = "Sample")
  expr_list <- lapply(names(my_study), function(sam) {
    message("  Processing sample: ", sam)
    seu <- my_study[[sam]]
    seu <- NormalizeData(seu, scale.factor = 1e5, verbose = FALSE)
    expr <- GetAssayData(seu, assay = "RNA", layer = "data") / log(2)
    expr
  })
  names(expr_list) <- names(my_study)

  mp_assignment <- assign_mp_labels(expr_list, mp_list, min_genes = min_genes / 50, min_score = min_score)
  mp_scores_per_sample <- mp_assignment$mp_scores
  assign_mp_per_cell <- mp_assignment$labels

  # Plot unfiltered barplot
  old_order <- c("un_assigned", "MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- c(setNames("un_assigned", "un_assigned"), setNames(paste0("MP_", 1:9), old_order[-1]))

  df <- lapply(assign_mp_per_cell, function(x) data.frame(MP = x))
  for (i in 1:length(df)) df[[i]]$sample <- names(df)[i]
  df <- do.call(rbind, df)
  df$MP <- as.character(df$MP)
  df$MP <- unname(mp_map[df$MP])
  df$MP <- factor(df$MP, levels = c("un_assigned", paste0("MP_", 1:9)))

  col_unfiltered <- c("grey", "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF",
                      "#3C5488FF", "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF")
  plot_mp_sample_barplot(df, luad_epi, file.path(output_fig_dir, "05_MP_tumor_sample.pdf"),
                         col_unfiltered, include_unassigned = TRUE)

  # Filter unassigned cells
  assign_mp_per_cell_filtered <- lapply(assign_mp_per_cell, function(x) x[x != "un_assigned"])

  # Plot filtered barplot
  df_filter <- lapply(assign_mp_per_cell_filtered, function(x) data.frame(MP = x))
  for (i in 1:length(df_filter)) df_filter[[i]]$sample <- names(df_filter)[i]
  df_filter <- do.call(rbind, df_filter)
  old_order_filt <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map_filt <- setNames(paste0("MP_", 1:9), old_order_filt)
  df_filter$MP <- as.character(df_filter$MP)
  df_filter$MP <- unname(mp_map_filt[df_filter$MP])
  df_filter$MP <- factor(df_filter$MP, levels = paste0("MP_", 1:9))

  col_filtered <- c("#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
                    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF")
  plot_mp_sample_barplot(df_filter, luad_epi, file.path(output_fig_dir, "06_MP_tumor_filter_sample.pdf"),
                         col_filtered, include_unassigned = FALSE)

  # Combine scores and labels
  mp_scores_per_cell <- do.call(rbind, mp_scores_per_sample)

  df_all <- lapply(assign_mp_per_cell, function(x) data.frame(MP = x))
  for (i in 1:length(df_all)) df_all[[i]]$sample <- names(df_all)[i]
  df_all <- do.call(rbind, df_all)

  df_filter_all <- lapply(assign_mp_per_cell_filtered, function(x) data.frame(MP = x))
  for (i in 1:length(df_filter_all)) df_filter_all[[i]]$sample <- names(df_filter_all)[i]
  df_filter_all <- do.call(rbind, df_filter_all)

  df_all$MP_anno <- ifelse(df_all$MP == "un_assigned", "un_assigned",
    ifelse(!rownames(df_all) %in% rownames(df_filter_all), "filter_MP_cells", as.character(df_all$MP)))

  mp_label_score <- cbind(mp_scores_per_cell, df_all)
  rownames(mp_label_score) <- sub("^[^.]*\\.", "", rownames(mp_label_score))
  mp_label_score <- mp_label_score[colnames(luad_epi), ]
  mp_label_score$sample <- NULL

  luad_epi@meta.data <- cbind(luad_epi@meta.data, mp_label_score)

  saveRDS(luad_epi, file = file.path(output_data_dir, "03_LUAD_Epi_assigned_MP.rds"))
  save(mp_label_score, assign_mp_per_cell, assign_mp_per_cell_filtered,
       file = file.path(output_data_dir, "03_Assign_MP_per_cell_filtered.RData"))

  # Final filtered object with CNV scores
  luad_epi <- subset(luad_epi, MP != "un_assigned")
  cnv_score <- readRDS(file.path(RESULT_1_DIR, "results_data", "infercnv", "09_CNV_score.rds"))
  colnames(luad_epi) <- gsub("-", "_", colnames(luad_epi))
  colnames(luad_epi) <- gsub("_", "", colnames(luad_epi))
  cnv_score <- cnv_score[colnames(luad_epi), ]
  luad_epi$CNV_score <- cnv_score$CNV_score

  saveRDS(luad_epi, file = file.path(output_data_dir, "04_LUAD_Epi_assigned_MP_final.rds"))

  # Plot pie charts
  plot_mp_stage_pies(luad_epi@meta.data[, c("stage", "MP_anno")],
                     file.path(output_fig_dir, "07_MP_tumor_cells_stage_pie.pdf"))

  message("MP signature identification complete.")
  invisible(list(Cluster_list = cluster_list, MP_list = mp_list, labels = assign_mp_per_cell))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP signature identification: ", e$message)
      stop(e)
    }
  )
}
