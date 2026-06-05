#' Ligand-Receptor Pair Prognosis Analysis
# Script: 10_lr_prognosis
#'
#' Performs Kaplan-Meier survival analysis for ligand-receptor gene pairs
#' across multiple bulk RNA-seq cohorts (TCGA + GEO datasets).
#' Outputs both four-group and High-High vs Low-Low comparison plots.
#'
#' @return None. PDF figures are written to the result figure directory.
#' @author Script refactored from original analysis pipeline.
#' @export
run_lr_prognosis <- function() {
# Auto-locate config.R
config_path <- "config.R"
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(dirname(dirname(getwd()))), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- "~/LUAD-scMetaPrograms/config.R"
}
source(config_path)
  library(survival)
  library(survminer)

  load(file.path(BULK_DIR, "TCGA_GEO7.rdata"))

  title <- c("TCGA", "GSE11969", "GSE13213", "GSE26939", "GSE30219", "GSE31210", "GSE41271", "GSE72094")

  gene_inter <- c("MDK", "NCL")
  gene_inter <- c("SCGB3A2", "MARCO")
  gene_inter <- c("ANGPTL4", "CDH5")
  gene_inter <- c("ANGPTL4", "CDH11")
  gene_inter <- c("PLAU", "PLAUR")

  # Four-group comparison
  Gene_group <- list()
  Gene_exp <- list()

  for (i in seq_along(expr_list)) {
    GSE_expr <- expr_list[[title[i]]]
    if (all(gene_inter %in% rownames(GSE_expr))) {
      exp <- GSE_expr[rownames(GSE_expr) %in% gene_inter, ]
      Gene_exp[[title[i]]] <- exp
      for (j in seq_len(nrow(exp))) {
        exp[j, ] <- ifelse(exp[j, ] >= quantile(as.numeric(exp[j, ]))[[3]], "High", "Low")
      }
      Gene_group[[title[i]]] <- as.data.frame(t(exp))
    } else {
      cat(sprintf("Note: Dataset [%s] does not contain all target genes, skipped.\n", i))
    }
  }

  names(clin_list) <- names(expr_list)
  clin_list <- clin_list[names(Gene_group)]

  km_plot <- list()
  for (i in seq_along(clin_list)) {
    clin_list[[i]]$risk_group <- paste(
      Gene_group[[i]][[gene_inter[1]]],
      Gene_group[[i]][[gene_inter[2]]], sep = " - "
    )
    if (length(unique(clin_list[[i]]$risk_group)) == 4) {
      km_fit <- survfit(Surv(survival_time, vital_status) ~ risk_group, data = clin_list[[i]])
      km_plot[[paste0(title[[i]], "_", gene_inter[1], "_", gene_inter[2])]] <- ggsurvplot(
        km_fit, data = clin_list[[i]],
        title = paste0(title[[i]], "_", gene_inter[1], "_", gene_inter[2]),
        font.x = 14, font.y = 14,
        font.tickslab = 12, pval = TRUE, pval.size = 4.5, size = 1.2,
        palette = c("#E41A1C", "grey", "lightgrey", "#377EB8"),
        legend = c(0.8, 0.8),
        legend.labs = c("High - High", "High - Low", "Low - High", "Low - Low"),
        legend.title = "group",
        tables.theme = theme_cleantable()
      )
    }
  }

  pdf(
    file.path(RESULT_5_DIR, "results_figure", paste0(gene_inter[1], "_", gene_inter[2], "_prognosis.pdf")),
    width = 15, height = 6, onefile = FALSE
  )
  arrange_ggsurvplots(km_plot, ncol = 4, nrow = 2)
  dev.off()

  # High-High vs Low-Low comparison only
  Gene_group <- list()
  Gene_exp <- list()

  for (i in seq_along(expr_list)) {
    GSE_expr <- expr_list[[title[i]]]
    if (all(gene_inter %in% rownames(GSE_expr))) {
      exp <- GSE_expr[rownames(GSE_expr) %in% gene_inter, ]
      Gene_exp[[title[i]]] <- exp
      for (j in seq_len(nrow(exp))) {
        exp[j, ] <- ifelse(exp[j, ] >= quantile(as.numeric(exp[j, ]))[[3]], "High", "Low")
      }
      Gene_group[[title[i]]] <- as.data.frame(t(exp))
    } else {
      cat(sprintf("Note: Dataset [%s] does not contain all target genes, skipped.\n", i))
    }
  }

  names(clin_list) <- names(expr_list)
  clin_list <- clin_list[names(Gene_group)]

  km_plot <- list()
  for (i in seq_along(clin_list)) {
    clin_list[[i]]$risk_group <- paste(
      Gene_group[[i]][[gene_inter[1]]],
      Gene_group[[i]][[gene_inter[2]]], sep = " - "
    )
    clin_list[[i]] <- clin_list[[i]][clin_list[[i]]$risk_group %in% c("High - High", "Low - Low"), ]
    if (length(unique(clin_list[[i]]$risk_group)) == 2) {
      km_fit <- survfit(Surv(survival_time, vital_status) ~ risk_group, data = clin_list[[i]])
      km_plot[[paste0(title[[i]], "_", gene_inter[1], "_", gene_inter[2])]] <- ggsurvplot(
        km_fit, data = clin_list[[i]],
        title = paste0(title[[i]], "_", gene_inter[1], "_", gene_inter[2]),
        font.x = 14, font.y = 14,
        font.tickslab = 12, pval = TRUE, pval.size = 4.5, size = 1.2,
        palette = c("#E41A1C", "#377EB8"),
        legend = c(0.8, 0.8),
        legend.labs = c("High - High", "Low - Low"),
        legend.title = "group",
        tables.theme = theme_cleantable()
      )
    }
  }

  pdf(
    file.path(RESULT_5_DIR, "results_figure", paste0(gene_inter[1], "_", gene_inter[2], "_prognosis_High_Low.pdf")),
    width = 15, height = 6, onefile = FALSE
  )
  arrange_ggsurvplots(km_plot, ncol = 4, nrow = 2)
  dev.off()
}

if (!interactive()) {
  run_lr_prognosis()
}
