# ============================================================================
# Script: 06_mp_stage_association
# Module: 02_meta_programs
# Description:
#   Validate MP scores across bulk RNA-seq stage cohorts using GSVA/ssGSEA
#   and generate boxplots comparing early vs advanced stages.
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
  library(ggpubr)
  library(GSVA)
  library(lemon)
  library(patchwork)

# ---- Helper Functions ----

#' Load MP gene lists.
#' @param rdata_path Path to RData file.
#' @return Named list of MP gene vectors.
load_mp_list <- function(rdata_path) {
  message("Loading MP list from: ", rdata_path)
  load(rdata_path)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_list <- MP_list[old_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))
  mp_list
}

#' Load and score a stage validation dataset.
#' @param file_path Path to RData file with GSE_expr and clinical.
#' @param mp_list MP gene lists.
#' @return List with gsva_matrix and clinical data frame.
score_stage_dataset <- function(file_path, mp_list) {
  message("Processing: ", basename(file_path))
  load(file_path)
  param <- ssgseaParam(exprData = as.matrix(GSE_expr), geneSets = mp_list)
  gsva_matrix <- gsva(param)
  list(gsva = gsva_matrix, clinical = clinical)
}

#' Generate stage comparison boxplots for a dataset.
#' @param scored_data List from score_stage_dataset.
#' @param mp_list MP gene lists.
#' @param dataset_name Dataset identifier.
#' @param stage_levels Factor levels for stage.
#' @param comparisons List of comparison pairs.
#' @param fill_colors Named vector of fill colors.
#' @return List of ggplot objects.
make_stage_boxplots <- function(scored_data, mp_list, dataset_name,
                                stage_levels, comparisons, fill_colors) {
  boxplots <- list()
  for (i in 1:length(mp_list)) {
    exp <- data.frame(ssgsea = scored_data$gsva[i, ],
                      stage = scored_data$clinical$stage)
    exp$stage <- factor(exp$stage, levels = stage_levels)
    p <- ggboxplot(exp, "stage", "ssgsea",
                   title = names(mp_list)[i], fill = "stage", outlier.shape = NA) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.title.x = element_blank()) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test") +
      scale_fill_manual(values = fill_colors)
    boxplots[[paste(names(mp_list)[i], dataset_name)]] <- p
  }
  boxplots
}

#' Generate TCGA stage comparison boxplots.
#' @param scored_data List with gsva and clinical.
#' @param mp_list MP gene lists.
#' @return List of ggplot objects.
make_tcga_boxplots <- function(scored_data, mp_list) {
  boxplots <- list()
  for (i in 1:length(mp_list)) {
    exp <- data.frame(ssgsea = scored_data$gsva[i, ],
                      stage = scored_data$clinical$stage)
    exp$stage <- factor(exp$stage, levels = c("normal", "early", "advance"))
    p <- ggboxplot(exp, "stage", "ssgsea", fill = "stage",
                   outlier.shape = NA, title = names(mp_list)[i]) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.title.x = element_blank()) +
      stat_compare_means(
        comparisons = list(c("early", "normal"), c("advance", "early"), c("advance", "normal")),
        method = "wilcox.test"
      ) +
      scale_fill_manual(values = c("normal" = "#fce38a", "early" = "#1F77B4", "advance" = "#D62728"))
    boxplots[[names(mp_list)[i]]] <- p
  }
  boxplots
}

#' Generate TCGA detailed stage boxplots.
#' @param scored_data List with gsva and clinical.
#' @param mp_list MP gene lists.
#' @return List of ggplot objects.
make_tcga_split_boxplots <- function(scored_data, mp_list) {
  clinical <- scored_data$clinical
  clinical$Stage <- case_when(
    clinical$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "I",
    clinical$ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "II",
    clinical$ajcc_pathologic_stage %in% c("Stage IIIA", "Stage IIIB") ~ "III",
    clinical$ajcc_pathologic_stage %in% c("Stage IV") ~ "IV"
  )
  clinical$Stage[which(clinical$stage == "normal")] <- "normal"
  clinical$Stage <- factor(clinical$Stage, levels = c("normal", "I", "II", "III", "IV"))

  group <- levels(factor(clinical$Stage))
  comp <- combn(group, 2)
  my_comparisons <- list()
  for (j in 1:ncol(comp)) {
    my_comparisons[[j]] <- comp[, j]
  }

  boxplots <- list()
  for (i in 1:length(mp_list)) {
    exp <- data.frame(ssgsea = scored_data$gsva[i, ], stage = clinical$Stage)
    p <- ggboxplot(exp, "stage", "ssgsea", fill = "stage",
                   outlier.shape = NA, title = names(mp_list)[i]) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.title.x = element_blank()) +
      stat_compare_means(method = "t.test", hide.ns = FALSE,
                         comparisons = my_comparisons, label = "p.adj.format",
                         p.adjust.method = "BH", vjust = 0.02, bracket.size = 0.6) +
      scale_fill_manual(values = c("#79B99D", "#4DBBD5FF", "pink", "#F39B7FFF", "#E64B35FF"))
    boxplots[[names(mp_list)[i]]] <- p
  }
  boxplots
}

# ---- Main Function ----

#' Run MP stage association analysis.
#' @param rdata_path Path to MP RData file.
#' @param output_fig_dir Directory for output figures.
#' @return List of all generated boxplots.
run_analysis <- function(
    rdata_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_list <- load_mp_list(rdata_path)

  stage_files <- c(
    file.path(BULK_DIR, "stage_vaild_data", "GSE13213_stage_data.Rdata"),
    file.path(BULK_DIR, "stage_vaild_data", "GSE26939_stage_data.Rdata"),
    file.path(BULK_DIR, "stage_vaild_data", "GSE41271_stage_data.Rdata"),
    file.path(BULK_DIR, "stage_vaild_data", "GSE72094_stage_data.Rdata"),
    file.path(BULK_DIR, "stage_vaild_data", "GSE11969_stage_data.Rdata")
  )
  dataset_names <- c("GSE13213", "GSE26939", "GSE41271", "GSE72094", "GSE11969")

  fill_colors <- c("early" = "#1F77B4", "advance" = "#D62728")
  comparisons <- list(c("advance", "early"))

  combined_list <- list()
  for (j in 1:length(stage_files)) {
    scored <- score_stage_dataset(stage_files[j], mp_list)
    bps <- make_stage_boxplots(scored, mp_list, dataset_names[j],
                               c("early", "advance"), comparisons, fill_colors)
    combined_list <- append(combined_list, bps)
  }

  pdf(file = file.path(output_fig_dir, "11_stage_vaild.pdf"), width = 30, height = 20, onefile = FALSE)
  print(wrap_plots(combined_list, ncol = 9, nrow = 5) +
    plot_layout(guides = "collect") & theme(legend.position = "top"))
  dev.off()
  message("Saved stage validation boxplots.")

  # TCGA
  tcga_file <- file.path(BULK_DIR, "TCGA", "TCGA_sur_data_stage.Rdata")
  tcga_scored <- score_stage_dataset(tcga_file, mp_list)
  tcga_boxplots <- make_tcga_boxplots(tcga_scored, mp_list)

  pdf(file = file.path(output_fig_dir, "12_MP_stage_TCGA.pdf"), width = 12, height = 12, onefile = FALSE)
  print(wrap_plots(tcga_boxplots, ncol = 3, nrow = 3) +
    plot_layout(guides = "collect") & theme(legend.position = "top"))
  dev.off()
  message("Saved TCGA stage boxplots.")

  tcga_split <- make_tcga_split_boxplots(tcga_scored, mp_list)
  pdf(file = file.path(output_fig_dir, "12_MP_stage_TCGA_split.pdf"), width = 15, height = 12, onefile = FALSE)
  print(wrap_plots(tcga_split, ncol = 3, nrow = 3) +
    plot_layout(guides = "collect") & theme(legend.position = "top"))
  dev.off()
  message("Saved TCGA split stage boxplots.")

  message("MP stage association analysis complete.")
  invisible(list(stage = combined_list, tcga = tcga_boxplots, tcga_split = tcga_split))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP stage association: ", e$message)
      stop(e)
    }
  )
}
