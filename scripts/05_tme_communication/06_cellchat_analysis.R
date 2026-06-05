# =============================================================================
# Script: 06_cellchat_analysis
# CellChat Analysis: Tumor-TME Communication
# =============================================================================
# This script loads LUAD tumor and TME Seurat objects, merges them, runs
# CellChat to infer intercellular communication, and generates visualizations
# for tumor sending/receiving signals.
# =============================================================================
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
library(CellChat)
library(Seurat)
library(future)
library(ggplot2)
library(reshape2)
library(patchwork)

options(future.globals.maxSize = 60 * 1024^3)
options(stringsAsFactors = FALSE)

# ---- Cell type ordering vectors (shared across functions) ----

TME_SUBTYPES <- c(
  "CD4+ Exhaustion", "CD4+ Treg", "CD8+ Exhaustion",
  "FCGR3A- NK", "FCGR3A+ NK", "IFIT3+ B", "H3+ B",
  "S100A12+ Mono", "SPP1+ Mac", "SELENOP+ Mac",
  "CXCL10+ Mac", "VEGFA+ Mac", "C3+ Mac", "LILRB2+ Mono",
  "CD1C+ DC", "LAMP3+ DC", "VWA1+ Endo", "IL7R+ Endo",
  "CCL21+ Endo", "COL10A1+ Fibro", "IGF1+ Fibro"
)

TUMOR_MP_LEVELS <- c("MP_5", "MP_6", "MP_3", "MP_4", "MP_7", "MP_1", "MP_2", "MP_8", "MP_9")

#' Load and prepare tumor epithelial cells
#'
#' Loads the LUAD epithelial object with meta-program (MP) assignments and
#' reorders MP factors to the desired level order.
#'
#' @param path Path to the RDS file containing LUAD_Epi_assigned_MP.
#' @return A Seurat object with MP reordered and stored in cell_subtype.
load_tumor_cells <- function(path) {
  tumor <- readRDS(path)
  tumor$MP <- as.character(tumor$MP)

  old_order <- c("MP_1", "MP_8", "MP_2", "MP_7", "MP_4", "MP_6", "MP_5", "MP_3", "MP_9")
  mp_map <- setNames(paste0("MP_", seq_along(old_order)), old_order)
  tumor$MP <- unname(mp_map[tumor$MP])
  tumor$MP <- factor(tumor$MP, levels = paste0("MP_", 1:9))

  tumor$cell_subtype <- tumor$MP
  return(tumor)
}

#' Load TME cell annotations
#'
#' Loads individual TME Seurat objects, subsets to the desired subtypes,
#' and stores the subtype identity in cell_subtype.
#'
#' @param t_nk_path Path to T_NK annotation RDS.
#' @param b_path Path to B cell annotation RDS.
#' @param myeloid_path Path to myeloid annotation RDS.
#' @param endo_path Path to endothelial annotation RDS.
#' @param fibro_path Path to fibroblast annotation RDS.
#' @return A named list of Seurat objects.
load_tme_cells <- function(t_nk_path, b_path, myeloid_path, endo_path, fibro_path) {
  t_nk <- readRDS(t_nk_path)
  t_nk <- subset(t_nk, T_NK_subtype %in% c(
    "FCGR3A+ NK", "FCGR3A- NK", "CD8+ Exhaustion",
    "CD4+ Exhaustion", "CD4+ Treg"
  ))
  t_nk$cell_subtype <- t_nk$T_NK_subtype

  b <- readRDS(b_path)
  b <- subset(b, B_subtype_subtype %in% c("H3+ B", "IFIT3+ B"))
  b$cell_subtype <- b$B_subtype_subtype

  myeloid <- readRDS(myeloid_path)
  myeloid <- subset(myeloid, Myeloid_subtype_subtype %in% c(
    "SELENOP+ Mac", "VEGFA+ Mac", "SPP1+ Mac",
    "CXCL10+ Mac", "C3+ Mac", "S100A12+ Mono",
    "LILRB2+ Mono", "CD1C+ DC", "LAMP3+ DC"
  ))
  myeloid$cell_subtype <- myeloid$Myeloid_subtype_subtype

  endo <- readRDS(endo_path)
  endo <- subset(endo, Endo_subtype_subtype %in% c(
    "IL7R+ Endo", "VWA1+ Endo", "CCL21+ Endo"
  ))
  endo$cell_subtype <- endo$Endo_subtype_subtype

  fibro <- readRDS(fibro_path)
  fibro <- subset(fibro, Fibro_subtype_subtype %in% c(
    "COL10A1+ Fibro", "IGF1+ Fibro"
  ))
  fibro$cell_subtype <- fibro$Fibro_subtype_subtype

  list(
    T_NK_anno = t_nk,
    B_anno = b,
    Myeloid_anno = myeloid,
    Endo_anno = endo,
    Fibro_anno = fibro
  )
}

#' Merge Seurat objects
#'
#' Merges tumor and TME Seurat objects into a single object and sets
#' the default identity to cell_subtype.
#'
#' @param tumor Tumor Seurat object.
#' @param tme_list Named list of TME Seurat objects.
#' @return A merged Seurat object.
merge_all_cells <- function(tumor, tme_list) {
  obj_list <- c(list(tumor), tme_list)
  merged <- merge(obj_list[[1]], obj_list[2:length(obj_list)]) %>% JoinLayers()
  Idents(merged) <- merged$cell_subtype
  return(merged)
}

#' Run CellChat analysis workflow
#'
#' Creates a CellChat object from the merged data, subsets the ligand-receptor
#' database, identifies over-expressed genes and interactions, computes
#' communication probabilities, and aggregates the network.
#'
#' @param merged_obj A merged Seurat object.
#' @param db_search Character string for CellChatDB subset (default "Secreted Signaling").
#' @return A CellChat object with computed results.
run_cellchat_workflow <- function(merged_obj, db_search = "Secreted Signaling") {
  cellchat <- createCellChat(object = merged_obj)

  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = db_search)
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}

#' Plot aggregate CellChat network
#'
#' Generates circle plots showing the number and weight of interactions
#' and saves them to a PDF.
#'
#' @param cellchat A CellChat object with aggregated net results.
#' @param out_dir Directory to save PDF outputs.
plot_aggregate_network <- function(cellchat, out_dir) {
  pdf(file.path(out_dir, "cellchat_aggregate_network.pdf"), width = 8, height = 6)

  par(mfrow = c(1, 1), xpd = TRUE)
  netVisual_circle(cellchat@net$count,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight,
                   vertex.weight = 20,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = "Interaction weights/strength")

  dev.off()
}

#' Extract tumor-TME communication matrices
#'
#' Extracts and melts send/receive weight matrices between tumor MPs
#' and TME cell subtypes. Returns formatted data frames for plotting.
#'
#' @param cellchat A CellChat object.
#' @return A list containing tsnw, trnw, tall, TME, and Tumor vectors.
extract_communication_matrices <- function(cellchat) {
  tumor_send_net_weight <- cellchat@net$weight[1:9, 10:30]
  tumor_recevier_net_weight <- cellchat@net$weight[10:30, 1:9]

  tsnw <- melt(tumor_send_net_weight)
  trnw <- melt(tumor_recevier_net_weight)
  tall <- rbind(tsnw, trnw)

  tsnw$Var1 <- factor(tsnw$Var1, levels = TUMOR_MP_LEVELS)
  tsnw$Var2 <- factor(tsnw$Var2, levels = rev(TME_SUBTYPES))
  trnw$Var2 <- factor(trnw$Var2, levels = TUMOR_MP_LEVELS)
  trnw$Var1 <- factor(trnw$Var1, levels = rev(TME_SUBTYPES))

  list(tsnw = tsnw, trnw = trnw, tall = tall)
}

#' Plot tumor-TME communication dotplots
#'
#' Generates side-by-side dotplots for tumor sending and receiving signals.
#'
#' @param tsnw Melted send matrix data frame.
#' @param trnw Melted receive matrix data frame.
#' @param tall Combined melted data frame for shared scales.
#' @param out_file Output PDF file path.
plot_tumor_tme_dotplots <- function(tsnw, trnw, tall, out_file) {
  g1 <- ggplot(tsnw, aes(x = Var1, y = Var2, color = value, size = value)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(position = "top") +
    geom_vline(xintercept = seq(1.5, length(unique(tsnw$Var1)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(tsnw$Var2)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
      breaks = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))

  g2 <- ggplot(trnw, aes(x = Var2, y = Var1, color = value, size = value)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(position = "top") +
    geom_vline(xintercept = seq(1.5, length(unique(trnw$Var2)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(trnw$Var1)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
      breaks = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.7, title = "Weight"),
           size = guide_legend("Weight"))

  G1 <- g1 + scale_radius(
    range = c(2, 6),
    limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE))
  ) + theme(legend.position = "none")

  G2 <- g2 + scale_radius(
    range = c(2, 6),
    limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
    labels = c("weight > 0.0", "weight > 0.2", "weight > 0.4", "weight > 0.6")
  )

  pdf(out_file, width = 8, height = 6)
  print(G1 + G2)
  dev.off()
}

#' Plot tumor-sending ligand-receptor interactions
#'
#' Generates a dotplot of ligand-receptor interactions sent from tumor MPs
#' to TME cells.
#'
#' @param cellchat A CellChat object.
#' @param out_file Output PDF file path.
plot_tumor_sending_lr <- function(cellchat, out_file) {
  cells.level <- levels(cellchat@idents)

  df.net <- subsetCommunication(cellchat, slot.name = "net",
                                sources.use = cells.level[1:9],
                                targets.use = cells.level[10:30],
                                thresh = 0.05)

  df.net$source <- paste0(df.net$source, " ", "state")
  df.net$prob.original <- df.net$prob
  df.net$prob <- -1 / log(df.net$prob)
  df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")

  tumor_state <- paste0(TUMOR_MP_LEVELS, " state")
  source_target <- paste(rep(tumor_state, times = length(TME_SUBTYPES)),
                         rep(TME_SUBTYPES, each = length(TUMOR_MP_LEVELS)),
                         sep = " -> ")
  df.net$source.target <- factor(df.net$source.target, levels = source_target)

  LR_sort <- unique(df.net$interaction_name_2)
  df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = rev(LR_sort))

  g <- ggplot(df.net, aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(df.net$source.target)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(df.net$interaction_name_2)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    theme(text = element_text(size = 10), plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6)) +
    ggtitle("All Communications") +
    theme(plot.title = element_text(hjust = 0.5))

  print(g)
  ggsave(out_file, width = 25, height = 10)
}

#' Plot tumor-receiving ligand-receptor interactions
#'
#' Generates a dotplot of ligand-receptor interactions received by tumor MPs
#' from TME cells.
#'
#' @param cellchat A CellChat object.
#' @param out_file Output PDF file path.
plot_tumor_receiving_lr <- function(cellchat, out_file) {
  cells.level <- levels(cellchat@idents)

  df.net <- subsetCommunication(cellchat, slot.name = "net",
                                sources.use = cells.level[10:30],
                                targets.use = cells.level[1:9],
                                thresh = 0.05)

  df.net$target <- paste0(df.net$target, " ", "state")
  df.net$prob.original <- df.net$prob
  df.net$prob <- -1 / log(df.net$prob)
  df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")

  tumor_state <- paste0(TUMOR_MP_LEVELS, " state")
  source_target <- paste(rep(TME_SUBTYPES, each = length(TUMOR_MP_LEVELS)),
                         rep(tumor_state, times = length(TME_SUBTYPES)),
                         sep = " -> ")
  df.net$source.target <- factor(df.net$source.target, levels = source_target)

  g <- ggplot(df.net, aes(x = source.target, y = interaction_name_2, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(df.net$source.target)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(df.net$interaction_name_2)) - 0.5, 1),
               lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))

  print(g)
  ggsave(out_file, width = 20, height = 10)
}

#' Main execution function
#'
#' Orchestrates the full CellChat analysis pipeline from data loading through
#' visualization.
main <- function() {
  tumor_path <- file.path(RESULT_2_DIR, "results_data", "04_LUAD_Epi_assigned_MP_final.rds")
  t_nk_path <- file.path(RESULT_5_DIR, "results_data", "02_T_NK_anno_final.rds")
  b_path <- file.path(RESULT_5_DIR, "results_data", "04_B_subtype_anno_final.rds")
  myeloid_path <- file.path(RESULT_5_DIR, "results_data", "06_Myeloid_subtype_anno_final.rds")
  endo_path <- file.path(RESULT_5_DIR, "results_data", "08_Endo_subtype_anno_final.rds")
  fibro_path <- file.path(RESULT_5_DIR, "results_data", "10_Fibro_subtype_anno_final.rds")

  out_data_dir <- file.path(RESULT_5_DIR, "results_data")
  out_fig_dir <- file.path(RESULT_5_DIR, "results_figure")

  if (!dir.exists(out_data_dir)) dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_fig_dir)) dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading tumor cells...")
  tumor <- load_tumor_cells(tumor_path)

  message("Loading TME cells...")
  tme_list <- load_tme_cells(t_nk_path, b_path, myeloid_path, endo_path, fibro_path)

  message("Merging objects...")
  merged <- merge_all_cells(tumor, tme_list)

  message("Running CellChat workflow...")
  cellchat <- run_cellchat_workflow(merged)

  message("Plotting aggregate network...")
  plot_aggregate_network(cellchat, out_fig_dir)

  message("Saving CellChat object...")
  saveRDS(cellchat, file.path(out_data_dir, "cellchat.rds"))

  message("Extracting communication matrices...")
  comm <- extract_communication_matrices(cellchat)

  message("Plotting tumor-TME dotplots...")
  plot_tumor_tme_dotplots(
    comm$tsnw, comm$trnw, comm$tall,
    file.path(out_fig_dir, "26_MP_TME_send_recept.pdf")
  )

  message("Plotting tumor-sending LR interactions...")
  plot_tumor_sending_lr(cellchat, file.path(out_fig_dir, "all_send.pdf"))

  message("Plotting tumor-receiving LR interactions...")
  plot_tumor_receiving_lr(cellchat, file.path(out_fig_dir, "all_recept.pdf"))

  message("CellChat analysis complete.")
}

if (!interactive()) {
  main()
}
