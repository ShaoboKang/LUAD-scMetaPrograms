#' VPS72 Characterization Across Bulk and Single-Cell Data
# Script: 02_vps72_characterization
#'
#' @description
#' Analyzes VPS72 expression across TCGA stages, survival cohorts,
#' cell-cycle correlations, and single-cell cell types.
#'
#' @return Invisibly returns NULL. Writes figures to RESULT_4_DIR.
#'
#' @examples
#' \dontrun{
#' run_vps72_characterization()
#' }
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
suppressMessages({
  library(ggpubr)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(tinyarray)
  library(survminer)
  library(survival)
  library(GSVA)
  library(purrr)
  library(ggsci)
  library(Seurat)
})

#' VPS72 Expression by Stage
#'
#' @param exp Expression data frame with VPS72, stage, and Stage columns.
#' @param out_pdf Output PDF path.
#' @return Patchwork plot object.
plot_vps72_stage <- function(exp, out_pdf) {
  exp$stage <- factor(exp$stage, levels = c("normal", "early", "advance"))
  p1 <- ggboxplot(exp, "stage", "VPS72", fill = "stage", outlier.shape = NA) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.title.x = element_blank()) +
    stat_compare_means(
      comparisons = list(c("early", "normal"), c("advance", "early"), c("advance", "normal")),
      method = "wilcox.test"
    ) +
    scale_fill_manual(values = c("#fce38a", "#1F77B4", "#D62728"))

  exp$Stage[which(exp$stage == "normal")] <- "Normal"
  exp$Stage <- case_when(
    exp$Stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
    exp$Stage %in% c("Stage IIA", "Stage IIB") ~ "Stage II",
    exp$Stage %in% c("Stage IIIA", "Stage IIIB") ~ "Stage III",
    exp$Stage %in% c("Stage IV") ~ "Stage IV",
    TRUE ~ "Normal"
  )

  exp$Stage <- factor(exp$Stage, levels = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV"))
  group <- levels(factor(exp$Stage))
  comp <- combn(group, 2)
  my_comparisons <- list()
  for (j in seq_len(ncol(comp))) {
    my_comparisons[[j]] <- comp[, j]
  }

  p2 <- ggboxplot(exp, "Stage", "VPS72", fill = "Stage", outlier.shape = NA) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_lancet() +
    stat_compare_means(
      method = "t.test",
      hide.ns = FALSE,
      comparisons = my_comparisons,
      label = "p.adj.format",
      p.adjust.method = "BH",
      vjust = 0.02,
      bracket.size = 0.6
    )

  pdf(out_pdf)
  print(p2)
  dev.off()

  invisible(p1 + p2)
}

#' VPS72 Survival Analysis Across Cohorts
#'
#' @param sur_dir Directory containing survival validation data.
#' @param out_pdf Output PDF path.
plot_vps72_survival <- function(sur_dir, out_pdf) {
  datasets <- c(
    "TCGA", "GSE13213", "GSE26939", "GSE30219",
    "GSE31210", "GSE41271", "GSE72094", "GSE11969"
  )
  dir_paths <- file.path(sur_dir, paste0(datasets, "_sur.Rdata"))
  dir_paths[1] <- file.path(sur_dir, "TCGA_sur.Rdata")

  combined_list <- list()
  for (j in seq_along(dir_paths)) {
    load(file = dir_paths[j])
    if ("VPS72" %in% rownames(GSE_expr)) {
      surv <- data.frame(
        row.names = colnames(GSE_expr),
        sample = rownames(clinical),
        event = clinical$vital_status,
        time = clinical$survival_time
      )
      group <- ifelse(
        GSE_expr["VPS72", ] >= quantile(as.matrix(GSE_expr["VPS72", ]))[[3]],
        "High", "Low"
      )
      surv <- cbind(surv, VPS72 = group)
      fit <- survfit(Surv(time, event) ~ VPS72, data = surv)
      combined_list[[datasets[j]]] <- ggsurvplot(
        fit, data = surv,
        palette = c("#E41A1C", "#377EB8"),
        pval = TRUE, size = 1.2, legend = "none", pval.size = 5,
        title = paste0(datasets[j], "_VPS72"), legend.lab = c("High", "Low")
      )
    }
  }

  pdf(out_pdf)
  arrange_ggsurvplots(
    combined_list,
    print = TRUE,
    title = "Survival curves",
    ncol = 4, nrow = 2,
    risk.table.height = 0.25
  )
  dev.off()
}

#' VPS72 Correlation with Cell Cycle Signatures
#'
#' @param bulk_rdata Path to combined bulk RNA-seq RData.
#' @param out_pdf Output PDF path.
plot_vps72_cellcycle <- function(bulk_rdata, out_pdf) {
  load(bulk_rdata)

  geneset <- list(
    G2M = cc.genes$g2m.genes,
    G1S = cc.genes$s.genes
  )

  combined_list <- list()
  for (j in names(expr_list)) {
    GSE_expr <- expr_list[[j]]
    if ("VPS72" %in% rownames(GSE_expr)) {
      param <- ssgseaParam(
        exprData = as.matrix(GSE_expr),
        geneSets = geneset
      )
      gsva_matrix <- gsva(param)
      gsva_matrix <- as.data.frame(t(rbind(gsva_matrix, GSE_expr[c("VPS72", "UBE2C"), ])))
      target_columns <- setdiff(colnames(gsva_matrix), "VPS72")

      plot_list <- map(target_columns, ~ {
        ggscatter(gsva_matrix, x = "VPS72", y = .x,
                  add = "reg.line", conf.int = TRUE,
                  add.params = list(fill = "lightgray")) +
          stat_cor(
            method = "pearson",
            label.x = min(gsva_matrix$VPS72, na.rm = TRUE),
            label.y = min(gsva_matrix[[.x]], na.rm = TRUE),
            hjust = 0, vjust = 0
          ) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(title = paste0("VPS72 vs ", .x, " in ", j), x = "VPS72", y = .x) +
          theme(plot.title = element_text(size = 10))
      })
      combined_list <- append(combined_list, plot_list)
    }
  }

  pdf(out_pdf, width = 18, height = 10)
  print(wrap_plots(combined_list, ncol = 6) + plot_annotation(title = "Correlation: VPS72 vs Cellcycle"))
  dev.off()
}

#' VPS72 Expression in Single-Cell Data
#'
#' @param cell_anno_rds Path to cell annotation RDS.
#' @param epi_anno_rds Path to epithelial annotation RDS.
#' @param out_pdf Output PDF path.
plot_vps72_sc_expression <- function(cell_anno_rds, epi_anno_rds, out_pdf) {
  LUAD_cell_anno_final <- readRDS(cell_anno_rds)
  p1 <- VlnPlot(LUAD_cell_anno_final, "VPS72", pt.size = 0, group.by = "cell_type")

  LUAD_epi_anno_final <- readRDS(epi_anno_rds)
  p2 <- VlnPlot(LUAD_epi_anno_final, "VPS72", pt.size = 0, group.by = "Epi_subtype")

  LUAD_cell_anno_final$detailed_cell_type <- as.character(LUAD_cell_anno_final$cell_type)
  LUAD_cell_anno_final$detailed_cell_type[
    match(colnames(LUAD_epi_anno_final), colnames(LUAD_cell_anno_final))
  ] <- as.character(LUAD_epi_anno_final$Epi_subtype)

  LUAD_cell_anno_final <- subset(LUAD_cell_anno_final, detailed_cell_type != "Epithelial")
  LUAD_cell_anno_final$detailed_cell_type <- factor(
    LUAD_cell_anno_final$detailed_cell_type,
    levels = c(
      "Tumor", "Ciliated", "Club", "AT1", "AT2",
      "T/NK", "B", "Plasma", "Myeloid", "Mast", "Fibroblast", "Endothelial"
    )
  )

  custom_colors <- c(
    "#E64B35FF", "#3C5488FF", "#00A087FF", "#4DBBD5FF",
    "#F39B7FFF", "#D0AFC4", "#4991C1", "#89558D",
    "#AFC2D9", "#435B95", "#79B99D", "#F5A623"
  )

  pdf(out_pdf, width = 6, height = 4)
  print(
    VlnPlot(LUAD_cell_anno_final, c("VPS72"), pt.size = 0, group.by = "detailed_cell_type") +
      scale_color_manual(values = custom_colors) +
      theme(legend.position = "none")
  )
  dev.off()
}

#' Run Full VPS72 Characterization
#'
#' @return Invisibly returns NULL.
run_vps72_characterization <- function() {
  message("[02_vps72_characterization] Starting analysis...")

  stage_rdata <- file.path(BULK_DIR, "TCGA", "TCGA_sur_data_stage.Rdata")
  load(stage_rdata)
  exp <- data.frame(
    VPS72 = t(mRNA_exp["VPS72", ]),
    stage = clinical$stage,
    Stage = clinical$ajcc_pathologic_stage
  )

  plot_vps72_stage(
    exp,
    file.path(RESULT_4_DIR, "results_figure", "08_VPS72_exp_in_TCGA_stage.pdf")
  )

  plot_vps72_survival(
    file.path(BULK_DIR, "sur_vaild_data"),
    file.path(RESULT_4_DIR, "results_figure", "09_VPS72_survival.pdf")
  )

  plot_vps72_cellcycle(
    file.path(BULK_DIR, "TCGA_GEO7.rdata"),
    file.path(RESULT_4_DIR, "results_figure", "07_VPS72_cellcycle.pdf")
  )

  plot_vps72_sc_expression(
    file.path(RESULT_1_DIR, "results_data", "06_LUAD_cell_anno_final.rds"),
    file.path(RESULT_1_DIR, "results_data", "08_LUAD_epi_anno_final.rds"),
    file.path(RESULT_4_DIR, "results_figure", "10_VPS72_exp_in_sc_celltype.pdf")
  )

  message("[02_vps72_characterization] Completed successfully.")
  invisible(NULL)
}

if (!interactive()) {
  run_vps72_characterization()
}
