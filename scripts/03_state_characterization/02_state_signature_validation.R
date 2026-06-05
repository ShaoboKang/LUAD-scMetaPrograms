# ============================================================================
# Script: 02_state_signature_validation
# Module: 03_state_characterization
# Description: Validate tumor state signatures with module scoring and violin plots.
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
library(ggrepel)
library(readxl)
library(patchwork)
library(clusterProfiler)

# ---- Helper Functions ----

#' Load Seurat object and remap MP labels to standard order.
#' @param input_path Path to the Seurat RDS file.
#' @return Seurat object with remapped MP identities.
load_and_remap_seurat <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path)
  }

  message("[INFO] Loading Seurat object from: ", input_path)
  seurat_obj <- tryCatch({
    readRDS(input_path)
  }, error = function(e) {
    stop("Failed to read RDS file: ", conditionMessage(e))
  })

  seurat_obj$MP <- as.character(seurat_obj$MP)

  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  seurat_obj$MP <- unname(mp_map[seurat_obj$MP])
  seurat_obj$MP <- factor(seurat_obj$MP, levels = paste0("MP_", 1:9))

  return(seurat_obj)
}

#' Load MSigDB gene sets and build custom pathway list.
#' @param msigdb_dir Directory containing MSigDB GMT and RData files.
#' @return Named list of gene sets.
load_gene_sets <- function(msigdb_dir) {
  message("[INFO] Loading MSigDB gene sets from: ", msigdb_dir)

  hallmark_path <- file.path(msigdb_dir, "m_hallmark.RData")
  if (!file.exists(hallmark_path)) {
    stop("Hallmark file not found: ", hallmark_path)
  }
  load(hallmark_path)

  kegg_path <- file.path(msigdb_dir, "kegg_hsa.gmt")
  go_path <- file.path(msigdb_dir, "c5.go.v2024.1.Hs.symbols.gmt")
  h_path <- file.path(msigdb_dir, "h.all.v2024.1.Hs.symbols.gmt")

  msig <- rbind(
    read.gmt(kegg_path),
    read.gmt(go_path),
    read.gmt(h_path)
  )

  msig <- msig[-grep("GOCC_", msig$term), ]
  msig <- msig[-grep("GOMF_", msig$term), ]
  msig <- msig[-grep("hsa", msig$term), ]

  # Debug: show MHC CLASS II related terms
  print(unique(msig$term[grep("MHC_CLASS_II", msig$term)]))

  gene_sets <- list()

  gene_sets[["G2M"]] <- cc.genes$g2m.genes
  gene_sets[["G1S"]] <- cc.genes$s.genes
  gene_sets[["TNFA_SIGNALING_VIA_NFKB"]] <- m_hallmark[["TNFA_SIGNALING_VIA_NFKB"]]
  gene_sets[["MHC_CLASS_II"]] <- msig$gene[which(msig$term == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II")]
  gene_sets[["CYTOPLASMIC_TRANSLATION_INITIATION"]] <- msig$gene[which(msig$term == "GOBP_FORMATION_OF_CYTOPLASMIC_TRANSLATION_INITIATION_COMPLEX")]
  gene_sets[["OXIDATIVE_PHOSPHORYLATION"]] <- m_hallmark[["OXIDATIVE_PHOSPHORYLATION"]]
  gene_sets[["CELL_ADHESION"]] <- msig$gene[which(msig$term == "GOBP_CELL_ADHESION")]
  gene_sets[["EPITHELIAL_MESENCHYMAL_TRANSITION"]] <- m_hallmark[["EPITHELIAL_MESENCHYMAL_TRANSITION"]]
  gene_sets[["HYPOXIA"]] <- m_hallmark[["HYPOXIA"]]
  gene_sets[["GLYCOLYSIS"]] <- m_hallmark[["GLYCOLYSIS"]]

  return(gene_sets)
}

#' Normalize data and add module scores for gene sets.
#' @param seurat_obj Seurat object.
#' @param gene_sets Named list of gene sets.
#' @return Seurat object with module scores in meta.data.
score_pathways <- function(seurat_obj, gene_sets) {
  message("[INFO] Adding module scores...")

  seurat_obj <- AddModuleScore(seurat_obj, features = gene_sets)

  # Rename newly added score columns to gene set names
  score_cols <- (ncol(seurat_obj@meta.data) - length(gene_sets) + 1):ncol(seurat_obj@meta.data)
  colnames(seurat_obj@meta.data)[score_cols] <- names(gene_sets)

  return(seurat_obj)
}

#' Plot violin plots of pathway module scores.
#' @param seurat_obj Seurat object with module scores.
#' @param gene_sets Named list of gene sets.
#' @param output_path Output PDF path.
#' @return NULL (invisibly).
plot_pathway_scores <- function(seurat_obj, gene_sets, output_path) {
  message("[INFO] Plotting pathway scores to: ", output_path)

  colors <- c(
    "#DC0000FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF",
    "#8491B4FF", "#F39B7FFF", "#7E6148FF", "#91D1C2FF"
  )

  plot_list <- list()  # Store individual MP plots

  for (i in names(gene_sets)) {
    plot_list[[i]] <- VlnPlot(
      seurat_obj,
      features = i,
      pt.size = 0,
      group.by = "MP",
      cols = colors
    ) +
      geom_boxplot(width = 0.2, col = "black", fill = "white") +
      NoLegend()
  }

  p_combined <- wrap_plots(plot_list, nrow = 2, ncol = 5) +
    plot_annotation(
      title = "Pathway Scores in nine tumor cell state",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  pdf(file = output_path, width = 20, height = 9)
  print(p_combined)
  dev.off()
  message("[INFO] Pathway scores plot saved.")

  invisible(NULL)
}

# ---- Main Analysis Functions ----

#' Main analysis function for state signature validation.
#' @param input_path Path to input Seurat RDS file.
#' @param output_dir Directory for output files.
#' @param msigdb_dir Directory containing MSigDB gene sets.
#' @return List containing the output PDF path.
run_analysis <- function(input_path, output_dir, msigdb_dir) {
  # Ensure output subdirectories exist
  if (!dir.exists(file.path(output_dir, "results_figure"))) {
    dir.create(file.path(output_dir, "results_figure"), recursive = TRUE, showWarnings = FALSE)
  }

  # Step 1: Load gene sets
  gene_sets <- load_gene_sets(msigdb_dir)

  # Step 2: Load and remap Seurat object
  seurat_obj <- load_and_remap_seurat(input_path)

  # Step 3: Normalize data
  seurat_obj <- NormalizeData(seurat_obj)

  # Step 4: Add module scores
  seurat_obj <- score_pathways(seurat_obj, gene_sets)

  # Step 5: Plot pathway scores
  output_path <- file.path(output_dir, "results_figure", "02_tumor_state_pathway_scores.pdf")
  plot_pathway_scores(seurat_obj, gene_sets, output_path)

  message("[INFO] Analysis complete.")
  return(list(pathway_plot = output_path))
}

# ---- Entry Point ----

if (!interactive()) {
  message("[INFO] Running analysis pipeline...")
  run_analysis()
} else {
  message("[INFO] Script loaded. Call run_analysis() to execute.")
}
