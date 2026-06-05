# ============================================================================
# Script: 09_riskscore_drug_correlation
# Module: 08_clinical
# Description: Spearman correlation between risk score and predicted drug
#   sensitivity (IC50), volcano plot, and scatter plots for selected drugs.
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

# ---- Helper Functions ----

#' Compute Spearman correlation per drug
#'
#' @param sen Data frame with drug IC50 columns and RS column.
#' @param drug_cols Character vector of drug column names.
#' @return Data frame with correlation results.
compute_drug_correlations <- function(sen, drug_cols) {
  cor_results <- lapply(drug_cols, function(drug) {
    x <- sen[[drug]]
    y <- sen$RS
    keep <- complete.cases(x, y)
    x <- x[keep]
    y <- y[keep]

    if (length(x) < 3) {
      return(data.frame(
        Drug = drug,
        cor = NA,
        p_value = NA,
        n = length(x)
      ))
    }

    test <- suppressWarnings(cor.test(x, y, method = "spearman"))
    data.frame(
      Drug = drug,
      cor = unname(test$estimate),
      p_value = test$p.value,
      n = length(x)
    )
  })

  cor_df <- bind_rows(cor_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      logP = -log10(p_adj),
      Direction = case_when(
        cor > 0 ~ "Positive",
        cor < 0 ~ "Negative",
        TRUE ~ "NS"
      )
    )
  return(cor_df)
}

# ---- Main Function ----

#' Run drug sensitivity correlation analysis
#'
#' @param drug_pred_path Path to DrugPredictions CSV.
#' @param riskscore_path Path to Model_riskscore RDS.
#' @param output_fig_dir Directory for figures.
#' @return List with correlation data frame and plots.
run_analysis <- function(
    drug_pred_path = file.path(RESULT_8_DIR, "results_data",
                               "calcPhenotype_Output", "DrugPredictions.csv"),
    riskscore_path = file.path(RESULT_8_DIR, "results_data", "Model_riskscore.rds"),
    output_fig_dir = file.path(RESULT_8_DIR, "results_figure")
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(ggpubr)

  if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading drug predictions and risk scores...")
  sen <- read.csv(drug_pred_path, row.names = 1, header = TRUE, check.names = FALSE)
  model_riskscore <- readRDS(riskscore_path)

  sen$RS <- model_riskscore$TCGA$RS[match(rownames(sen), rownames(model_riskscore$TCGA))]
  message("Unmatched samples: ", sum(is.na(sen$RS)))

  drug_cols <- setdiff(colnames(sen), "RS")

  message("Computing drug-risk score correlations...")
  cor_df <- compute_drug_correlations(sen, drug_cols)

  label_df <- cor_df %>%
    filter(!is.na(p_adj), p_adj < 0.05) %>%
    mutate(rank_score = abs(cor) * logP) %>%
    arrange(desc(rank_score)) %>%
    slice(1:20)

  cor_df <- cor_df %>%
    mutate(
      Sig = case_when(
        p_adj < 0.05 & cor > 0.3 ~ "Positive",
        p_adj < 0.05 & cor < -0.3 ~ "Negative",
        TRUE ~ "NS"
      )
    )

  p_volcano <- ggplot(cor_df, aes(x = cor, y = logP)) +
    geom_point(aes(color = Sig), size = 2.2, alpha = 0.8) +
    scale_color_manual(
      values = c(
        "Positive" = "#D55E00",
        "Negative" = "#0072B2",
        "NS" = "grey75"
      )
    ) +
    geom_vline(xintercept = c(-0.3, 0, 0.3), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = label_df,
      aes(label = Drug),
      size = 3.5,
      segment.size = 0.5,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey40",
      max.overlaps = Inf
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 13) +
    labs(
      x = "Spearman correlation with riskscore",
      y = expression(-log[10]("p.adj")),
      color = NULL
    ) +
    theme(
      legend.position = "top",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(10, 20, 30, 20)
    )

  print(p_volcano)

  pdf(file.path(output_fig_dir, "DrugSensitivity_riskscore_volcano.pdf"),
      width = 8, height = 6)
  print(p_volcano)
  dev.off()

  # Correlation scatter plots for selected drugs
  message("Generating scatter plots for selected drugs...")
  sel_drugs <- c("nutlin-3", "nilotinib", "bexarotene",
                 "MK-2206", "MK-1775", "paclitaxel")
  missing_drugs <- setdiff(sel_drugs, colnames(sen))
  if (length(missing_drugs) > 0) {
    message("Missing drugs: ", paste(missing_drugs, collapse = ", "))
  }

  plot_df <- sen %>%
    rownames_to_column("Sample_ID") %>%
    select(Sample_ID, RS, all_of(sel_drugs)) %>%
    pivot_longer(
      cols = all_of(sel_drugs),
      names_to = "Drug",
      values_to = "Predicted_IC50"
    ) %>%
    filter(!is.na(RS), !is.na(Predicted_IC50))

  plot_df$Drug <- factor(plot_df$Drug, levels = sel_drugs)

  p_scatter <- ggplot(plot_df, aes(x = RS, y = Predicted_IC50)) +
    geom_point(size = 1.2, alpha = 0.7, color = "#4C78A8") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.8) +
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 4
    ) +
    facet_wrap(~ Drug, scales = "free_y", ncol = 6) +
    theme_classic(base_size = 13) +
    labs(x = "Risk score", y = "Predicted IC50") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  print(p_scatter)

  pdf(file.path(output_fig_dir, "Drug_RS_correlation_scatter.pdf"),
      width = 18, height = 3.5)
  print(p_scatter)
  dev.off()

  message("Drug sensitivity correlation analysis complete.")
  return(list(cor_df = cor_df, volcano = p_volcano, scatter = p_scatter))
}

# ---- Entry Point ----
if (!interactive()) {
  results <- run_analysis()
}
