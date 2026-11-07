# Manuscript Methods 2.2.13; Results 3.6, Fig. S9.
# Adapted from bulk diagnostic Step_ssgsea_MP_AUC.R; no model fitting.
source("config.R")
suppressPackageStartupMessages({library(GSVA); library(pROC); library(ggplot2)})

#' Score all nine legacy MP signatures in seven diagnostic cohorts.
#' @param input_path Serialized list of sample-by-gene tables (TYPE Y/N).
#' @param output_dir Directory for scores, AUC tables and plots.
#' @return Long AUC table, invisibly.
run_analysis <- function(input_path = file.path(DIAGNOSTIC_DIR, "all_datasets_merged_list.rds"),
                         output_dir = file.path(RESULT_8_DIR, "diagnostic")) {
  e <- new.env()
  load(file.path(RESULT_2_DIR, "results_data", "02_Malig_final_MP_top50.RData"), e)
  stopifnot(all(MP_LEGACY_ORDER %in% names(e$MP_list)))
  sets <- lapply(e$MP_list[MP_LEGACY_ORDER], function(g) unique(toupper(g)))
  names(sets) <- MP_STATE_NAMES
  datasets <- readRDS(input_path)
  expected <- c("TCGA", "GSE31210", "GSE10072", "GSE19188", "GSE32863", "GSE63459", "GSE68571")
  stopifnot(setequal(names(datasets), expected))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  scores <- list(); curves <- list(); rows <- list()
  for (ds in expected) {
    d <- datasets[[ds]]
    names(d) <- toupper(trimws(names(d)))
    stopifnot("TYPE" %in% names(d), all(d$TYPE %in% c("Y", "N")))
    group <- factor(ifelse(d$TYPE == "Y", "Tumor", "Normal"), levels = c("Normal", "Tumor"))
    stopifnot(all(table(group) > 0))
    gene_cols <- which(!names(d) %in% c("SAMPLE", "TYPE", "STAGE") & !duplicated(names(d)))
    expr <- t(vapply(d[, gene_cols, drop = FALSE], function(x) as.numeric(as.character(x)), numeric(nrow(d))))
    stopifnot(all(is.finite(expr)))
    if (ds == "GSE68571") expr <- log2(pmax(expr, 0) + 1)
    coverage <- vapply(sets, function(g) sum(g %in% rownames(expr)), integer(1))
    stopifnot(all(coverage > 0))
    score <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr, geneSets = sets), verbose = FALSE)
    scores[[ds]] <- score
    for (mp in names(sets)) {
      # Fixed direction matches the original analysis; AUC below 0.5 is retained.
      r <- pROC::roc(group, score[mp, ], levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
      ci <- pROC::ci.auc(r)
      key <- paste(ds, mp, sep = "__")
      rows[[key]] <- data.frame(Dataset = ds, State = mp, AUC = as.numeric(pROC::auc(r)),
        CI_lower = ci[1], CI_upper = ci[3], Genes_used = coverage[mp],
        Normal = sum(group == "Normal"), Tumor = sum(group == "Tumor"))
      curves[[key]] <- data.frame(Dataset = ds, State = mp,
        FPR = 1 - r$specificities, TPR = r$sensitivities)
    }
  }
  auc <- do.call(rbind, rows); rownames(auc) <- NULL
  roc <- do.call(rbind, curves)
  saveRDS(scores, file.path(output_dir, "MP_ssgsea_scores_all_datasets.rds"))
  write.csv(auc, file.path(output_dir, "MP_diagnostic_AUC.csv"), row.names = FALSE)
  auc$State <- factor(auc$State, levels = rev(MP_STATE_NAMES))
  p <- ggplot(auc, aes(Dataset, State, fill = AUC)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", AUC))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = .5, limits = c(0, 1)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_dir, "MP_diagnostic_AUC_heatmap.pdf"), p, width = 9, height = 6)
  p <- ggplot(roc, aes(FPR, TPR, color = State)) + geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + facet_wrap(~Dataset) + theme_bw()
  ggsave(file.path(output_dir, "MP_diagnostic_ROC.pdf"), p, width = 12, height = 9)
  invisible(auc)
}
if (!interactive()) run_analysis()
