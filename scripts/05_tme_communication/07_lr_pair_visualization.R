#' Ligand-Receptor Pair Visualization
# Script: 07_lr_pair_visualization
#'
#' Generates dot plots, box plots, bar plots, and violin plots for selected
#' ligand-receptor pairs using CellChat and Seurat objects.
#'
#' @return None. PDF figures are written to the result figure directory.
#' @author Script refactored from original analysis pipeline.
#' @export
run_lr_pair_visualization <- function() {
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
  library(patchwork)
  options(stringsAsFactors = FALSE)
  library(mindr)
  library(Seurat)
  library(ggalluvial)

  cellchat <- readRDS(file.path(RESULT_5_DIR, "results_data", "cellchat.rds"))
  cells.level <- levels(cellchat@idents)

  # ggplot2 theme settings
  tttheme <- theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  ) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, vjust = 0.5)
    )

  # Define cell type order
  TME <- c(
    "CD4+ Exhaustion", "CD4+ Treg", "CD8+ Exhaustion", "FCGR3A- NK", "FCGR3A+ NK",
    "IFIT3+ B", "H3+ B", "S100A12+ Mono", "SPP1+ Mac", "SELENOP+ Mac",
    "CXCL10+ Mac", "VEGFA+ Mac", "C3+ Mac", "LILRB2+ Mono", "CD1C+ DC",
    "LAMP3+ DC", "VWA1+ Endo", "IL7R+ Endo", "CCL21+ Endo", "COL10A1+ Fibro",
    "IGF1+ Fibro"
  )

  Tumor <- c("MP_5", "MP_6", "MP_3", "MP_4", "MP_7", "MP_1", "MP_2", "MP_8", "MP_9")

  # Violin plot data preparation
  w10x <- CreateSeuratObject(cellchat@data.signaling, meta.data = cellchat@meta)
  Idents(w10x) <- w10x$cell_subtype
  w10x_tumor <- subset(w10x, idents = Tumor)
  Idents(w10x_tumor) <- factor(w10x_tumor$cell_subtype, levels = Tumor)
  w10x_TME <- subset(w10x, idents = TME)
  Idents(w10x_TME) <- factor(w10x_TME$cell_subtype, levels = TME)

  tumor_color <- c(
    "#3C5488FF", "#8491B4FF", "#00A087FF", "#4DBBD5FF",
    "#F39B7FFF", "#DC0000FF", "#E64B35FF", "#7E6148FF", "#91D1C2FF"
  )

  # Common LR pair MDK-NCL
  pair <- data.frame(interaction_name = "MDK_NCL")
  MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                            sources.use = cells.level[1:9], targets.use = cells.level[10:30],
                            thresh = 0.05)
  MN$prob <- -1 / log(MN$prob)
  MN$source <- factor(MN$source, levels = rev(Tumor))
  MN$target <- factor(MN$target, levels = TME)
  pair_all <- MN

  g1 <- ggplot(MN, aes(x = source, y = target, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank()) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(MN$source)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(MN$target)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(MN$prob, 0, na.rm = TRUE), quantile(MN$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(MN$prob, 0, na.rm = TRUE), quantile(MN$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    guides(size = guide_legend("Commun. Prob.")) +
    xlab("") + ylab("") +
    labs(title = "MDK - NCL") + tttheme + theme(legend.position = "none") +
    coord_flip()

  x <- MN
  x$source <- factor(x$source, levels = Tumor)
  g2 <- ggplot(x) +
    geom_boxplot(aes(source, prob, color = source)) +
    geom_jitter(aes(source, prob, color = source), width = 0.01) +
    scale_color_manual(values = tumor_color) +
    labs(title = "MDK - NCL") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    tttheme +
    theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

  g3 <- VlnPlot(w10x_tumor, features = "MDK", pt.size = 0, cols = tumor_color) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  g4 <- VlnPlot(w10x_TME, features = "NCL", pt.size = 0) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  layout <- "AB
CC"
  pdf(file.path(RESULT_5_DIR, "results_figure", "MDK_NCL_exp.pdf"), width = 8, height = 6)
  g2 + g3 + g4 + plot_layout(design = layout)
  dev.off()

  # Common LR pair SCGB3A2-MARCO
  pair <- data.frame(interaction_name = "SCGB3A2_MARCO")
  MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                            sources.use = cells.level[1:9], targets.use = cells.level[10:30],
                            thresh = 0.05)
  MN$prob <- -1 / log(MN$prob)
  MN$source <- factor(MN$source, levels = rev(Tumor))
  MN$target <- factor(MN$target, levels = TME)
  pair_all <- MN

  x <- MN
  x$source <- factor(x$source, levels = Tumor)
  g2 <- ggplot(x) +
    geom_boxplot(aes(source, prob, color = source)) +
    geom_jitter(aes(source, prob, color = source), width = 0.01) +
    scale_color_manual(values = tumor_color) +
    labs(title = "SCGB3A2 - MARCO") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    tttheme +
    theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

  g3 <- VlnPlot(w10x_tumor, features = "SCGB3A2", pt.size = 0, cols = tumor_color) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  g4 <- VlnPlot(w10x_TME, features = "MARCO", pt.size = 0) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  layout <- "ABCC"
  pdf(file.path(RESULT_5_DIR, "results_figure", "SCGB3A2_MARCO_exp.pdf"), width = 12, height = 4)
  g2 + g3 + g4 + plot_layout(design = layout)
  dev.off()

  # Early LR pair ANGPTL4-CDH5
  pair <- data.frame(interaction_name = "ANGPTL4_CDH5")
  MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                            sources.use = cells.level[1:9], targets.use = cells.level[10:30],
                            thresh = 0.05)
  MN$prob <- -1 / log(MN$prob)
  MN$source <- factor(MN$source, levels = rev(Tumor))
  MN$target <- factor(MN$target, levels = TME)

  x <- MN
  x$LR <- paste0(x$source, x$target)

  pair_all$LR <- paste0(pair_all$source, pair_all$target)
  pair_all_1 <- pair_all[-match(x$LR, pair_all$LR), ]
  pair_all_1$prob <- NA
  x <- rbind(x, pair_all_1)
  x$source <- factor(x$source, levels = Tumor)

  g1 <- ggplot(x, aes(x = source, y = target, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank()) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(x$source)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(x$target)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    guides(size = guide_legend("Commun. Prob.")) +
    xlab("") + ylab("") +
    labs(title = "ANGPTL4_CDH5") + tttheme + theme(legend.position = "none") + coord_flip()

  gdotplot1 <- g1

  x_mean <- x %>%
    group_by(source) %>%
    summarise(prob = mean(prob, na.rm = TRUE)) %>%
    ungroup()

  g2 <- ggplot(x_mean, aes(source, prob, fill = source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = tumor_color) +
    labs(title = "ANGPTL4_CDH5") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    tttheme +
    theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

  g3 <- VlnPlot(w10x_tumor, features = "ANGPTL4", pt.size = 0, cols = tumor_color) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")
  g4 <- VlnPlot(w10x_TME, features = "CDH5", pt.size = 0) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  layout <- "AAABBCCDDD"
  pdf(file.path(RESULT_5_DIR, "results_figure", "ANGPTL4_CDH5_exp.pdf"), width = 16, height = 4)
  g1 + g2 + g3 + g4 + plot_layout(design = layout)
  dev.off()

  # LR pair ANGPTL4-CDH11
  pair <- data.frame(interaction_name = "ANGPTL4_CDH11")
  MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                            sources.use = cells.level[1:9], targets.use = cells.level[10:30],
                            thresh = 0.05)
  MN$prob <- -1 / log(MN$prob)
  MN$source <- factor(MN$source, levels = rev(Tumor))
  MN$target <- factor(MN$target, levels = TME)
  x <- MN
  x$LR <- paste0(x$source, x$target)

  pair_all$LR <- paste0(pair_all$source, pair_all$target)
  pair_all_1 <- pair_all[-match(x$LR, pair_all$LR), ]
  pair_all_1$prob <- NA
  x <- rbind(x, pair_all_1)
  x$source <- factor(x$source, levels = Tumor)

  g1 <- ggplot(x, aes(x = source, y = target, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank()) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(x$source)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(x$target)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    guides(size = guide_legend("Commun. Prob.")) +
    xlab("") + ylab("") +
    labs(title = "ANGPTL4_CDH11") + tttheme + theme(legend.position = "none") + coord_flip()

  gdotplot3 <- g1

  x$source <- factor(x$source, levels = Tumor)
  x_mean <- x %>%
    group_by(source) %>%
    summarise(prob = mean(prob, na.rm = TRUE)) %>%
    ungroup()

  g2 <- ggplot(x_mean, aes(source, prob, fill = source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = tumor_color) +
    labs(title = "ANGPTL4_CDH11") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    tttheme +
    theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

  g3 <- VlnPlot(w10x_TME, features = c("CDH11"), pt.size = 0) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  layout <- "AAABBCCC"
  pdf(file.path(RESULT_5_DIR, "results_figure", "ANGPTL4_CDH511_exp.pdf"), width = 12, height = 4)
  g1 + g2 + g3 + plot_layout(design = layout)
  dev.off()

  # LR pair PLAU-PLAUR
  pair <- data.frame(interaction_name = "PLAU_PLAUR")
  MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                            sources.use = cells.level[1:9], targets.use = cells.level[10:30],
                            thresh = 0.05)
  MN$prob <- -1 / log(MN$prob)
  MN$source <- factor(MN$source, levels = rev(Tumor))
  MN$target <- factor(MN$target, levels = TME)
  x <- MN
  x$LR <- paste0(x$source, x$target)

  pair_all$LR <- paste0(pair_all$source, pair_all$target)
  pair_all_1 <- pair_all[-match(x$LR, pair_all$LR), ]
  pair_all_1$prob <- NA
  x <- rbind(x, pair_all_1)
  x$source <- factor(x$source, levels = Tumor)

  g5 <- ggplot(x, aes(x = source, y = target, color = prob, size = prob)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank()) +
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(1.5, length(unique(x$source)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    geom_hline(yintercept = seq(1.5, length(unique(x$target)) - 0.5, 1), lwd = 0.1, colour = "grey90") +
    scale_colour_gradientn(
      colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99),
      na.value = "white",
      limits = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      breaks = c(quantile(x$prob, 0, na.rm = TRUE), quantile(x$prob, 1, na.rm = TRUE)),
      labels = c("min", "max")
    ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    guides(size = guide_legend("Commun. Prob.")) +
    xlab("") + ylab("") +
    labs(title = "PLAU_PLAUR") + tttheme + theme(legend.position = "none") + coord_flip()

  x_mean <- x %>%
    group_by(source) %>%
    summarise(prob = mean(prob, na.rm = TRUE)) %>%
    ungroup()

  g6 <- ggplot(x_mean, aes(source, prob, fill = source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = tumor_color) +
    labs(title = "PLAU_PLAUR") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    tttheme +
    theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

  g7 <- VlnPlot(w10x_tumor, features = "PLAU", pt.size = 0, cols = tumor_color) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")
  g8 <- VlnPlot(w10x_TME, features = c("PLAUR"), pt.size = 0) +
    theme_classic() + tttheme + theme(legend.position = "none") + xlab("") +
    ylab("Expression")

  layout <- "AAABBCCDDD"
  pdf(file.path(RESULT_5_DIR, "results_figure", "PLAU_PLAUR_exp.pdf"), width = 16, height = 4)
  g5 + g6 + g7 + g8 + plot_layout(design = layout)
  dev.off()
}

if (!interactive()) {
  run_lr_pair_visualization()
}
