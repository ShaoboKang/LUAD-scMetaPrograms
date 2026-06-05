# ============================================================================
# Script: 04_mp_functional_annotation
# Module: 02_meta_programs
# Description:
#   Re-annotate MP genes by GO, KEGG, and Hallmark pathway enrichment,
#   and generate bubble plots of top enriched terms.
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
  library(clusterProfiler)
  library(readxl)
  library(viridis)
  library(parallel)
  library(ggplot2)
  library(scales)
  library(forcats)

# ---- Helper Functions ----

#' Load MP gene lists.
#' @param rdata_path Path to RData file with MP_list.
#' @return Named list of MP gene vectors.
load_mp_list <- function(rdata_path) {
  message("Loading MP list from: ", rdata_path)
  load(rdata_path)
  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_list <- MP_list[old_order]
  names(mp_list) <- paste0("MP_", 1:length(mp_list))
  mp_list
}

#' Load and combine MSigDB gene sets.
#' @param msigdb_dir Directory containing GMT files.
#' @return Combined data frame of gene sets.
load_msigdb_sets <- function(msigdb_dir) {
  message("Loading MSigDB gene sets...")
  msig <- rbind(
    read.gmt(file.path(msigdb_dir, "kegg_hsa.gmt")),
    read.gmt(file.path(msigdb_dir, "c5.go.v2024.1.Hs.symbols.gmt")),
    read.gmt(file.path(msigdb_dir, "h.all.v2024.1.Hs.symbols.gmt"))
  )
  msig <- msig[-grep("GOCC_", msig$term), ]
  msig <- msig[-grep("GOMF_", msig$term), ]
  msig <- msig[-grep("hsa", msig$term), ]
  msig
}

#' Run enrichment analysis for each MP.
#' @param mp_list Named list of gene vectors.
#' @param msig MSigDB data frame.
#' @param n_cores Number of cores. Default 10.
#' @return Named list of enrichment results.
enrich_mp_pathways <- function(mp_list, msig, n_cores = 10) {
  cl <- makeCluster(n_cores, type = "FORK")
  mp_msig <- parLapply(cl = cl, mp_list, function(x) {
    en <- enricher(x, TERM2GENE = msig, minGSSize = 5, maxGSSize = 2000)
    en@result
  })
  stopCluster(cl)
  mp_msig <- lapply(mp_msig, function(x) x[x$p.adjust < 0.05, ])
  mp_msig
}

#' Extract top pathways for plotting.
#' @param mp_msig List of enrichment results.
#' @param top_n Number of top terms per MP. Default 5.
#' @return Combined data frame of selected pathways.
extract_top_pathways <- function(mp_msig, top_n = 5) {
  pathway <- data.frame()
  for (i in 1:length(mp_msig)) {
    if (dim(mp_msig[[i]])[1] > top_n - 1) {
      tmp <- mp_msig[[i]][1:top_n, ]
    } else {
      tmp <- mp_msig[[i]]
    }
    tmp$MP <- names(mp_msig)[i]
    pathway <- rbind(pathway, tmp)
  }
  pathway
}

#' Shorten long pathway names.
#' @param pathway Data frame with Description column.
#' @return Data frame with shortened descriptions.
shorten_pathway_names <- function(pathway) {
  pathway$Description[which(pathway$Description == "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX")] <- "GOBP_PEPTIDE_MHC_COMPLEX_ASSEMBLY"
  pathway$Description[which(pathway$Description == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN")] <- "GOBP_EXOG_PEPTIDE_AG_PROC_PRES"
  pathway$Description[which(pathway$Description == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II")] <- "GOBP_MHCII_PEPTIDE_OR_POLYSACCH_AG_PROC_PRES"
  pathway$Description[which(pathway$Description == "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX")] <- "GOBP_PEPTIDE_MHCII_COMPLEX_ASSEMBLY"
  pathway$Description[which(pathway$Description == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II")] <- "GOBP_MHCII_EXOG_PEPTIDE_AG_PROC_PRES"
  pathway
}

#' Plot enrichment bubble chart.
#' @param pathway Data frame with MP, Description, p.adjust, Count.
#' @param output_path Output PDF path.
#' @param width PDF width.
#' @param height PDF height.
plot_enrichment_bubble <- function(pathway, output_path, width = 7, height = 7) {
  pathway$Description <- as.factor(pathway$Description)
  pathway$Description <- fct_inorder(pathway$Description)
  pathway$GeneRatio <- sapply(pathway$GeneRatio, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    parts[1] / parts[2]
  })
  pathway$MP <- factor(pathway$MP, levels = c(paste0("MP_", 1:9)))

  p <- ggplot(pathway, aes(MP, Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_viridis_c(option = "turbo", direction = -1, name = "p.adjust") +
    labs(x = NULL, y = NULL) +
    guides(size = guide_legend(order = 1)) +
    theme(legend.direction = "horizontal", legend.position = "bottom") +
    scale_y_discrete(position = "right", limits = rev(levels(pathway$Description)))

  pdf(file = output_path, width = width, height = height)
  print(p)
  dev.off()
  message("Saved enrichment bubble plot to: ", output_path)
}

# ---- Main Function ----

#' Run MP functional annotation.
#' @param rdata_path Path to MP RData file.
#' @param msigdb_dir Directory with MSigDB GMT files.
#' @param output_data_dir Directory for output data.
#' @param output_fig_dir Directory for output figures.
#' @return List of enrichment results.
run_analysis <- function(
    rdata_path = file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"),
    msigdb_dir = MSIGDB_DIR,
    output_data_dir = file.path(RESULT_2_DIR, "results_data"),
    output_fig_dir = file.path(RESULT_2_DIR, "results_figure")
) {
  if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_list <- load_mp_list(rdata_path)
  msig <- load_msigdb_sets(msigdb_dir)

  message("Running pathway enrichment...")
  mp_msig <- enrich_mp_pathways(mp_list, msig, n_cores = 10)

  pathway <- extract_top_pathways(mp_msig, top_n = 5)

  # Move HALLMARK_HYPOXIA and HALLMARK_GLYCOLYSIS to end
  hyp_idx <- which(pathway$Description == "HALLMARK_HYPOXIA")
  gly_idx <- which(pathway$Description == "HALLMARK_GLYCOLYSIS")
  other_idx <- setdiff(1:nrow(pathway), c(hyp_idx, gly_idx))
  pathway <- pathway[c(other_idx, hyp_idx, gly_idx), ]

  pathway <- shorten_pathway_names(pathway)

  plot_enrichment_bubble(pathway, file.path(output_fig_dir, "08_MP_function_define_go.pdf"))
  plot_enrichment_bubble(pathway, file.path(output_fig_dir, "08_MP_function_define_go_custom_colors.pdf"))

  # write.csv(pathway[, -c(1:2)], file = file.path(output_data_dir, "06_MP_pathway.csv"))

  message("MP functional annotation complete.")
  invisible(mp_msig)
}

# ---- Entry Point ----
if (!interactive()) {
  results <- tryCatch(
    run_analysis(),
    error = function(e) {
      message("Error in MP functional annotation: ", e$message)
      stop(e)
    }
  )
}
