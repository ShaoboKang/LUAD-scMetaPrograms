# Methods 2.2.15; Fig. 7A and S10A. Exact interventional SHAP for Cox log risk.
# Adapted from SHAP.R: fixed model, original inputs, fixed TCGA background.
source("config.R")
suppressPackageStartupMessages({library(glmnet); library(ggplot2)})

#' Explain the saved eight-TF model without refitting or changing preprocessing.
#' @param model_dir Folder containing Mime result and exact post-ComBat inputs.
#' @param output_dir Directory for audit tables and SHAP plots.
#' @return List of SHAP matrices and per-cohort numerical audit.
run_analysis <- function(model_dir = MODEL_DIR, output_dir = SHAP_DIR) {
  res <- readRDS(file.path(model_dir, "Mime_res_batch_corrected.rds"))
  inputs <- readRDS(file.path(model_dir, "data_list_for_SHAP.rds"))
  fit <- res$ml.res[[FINAL_MODEL]]
  scores <- res$riskscore[[FINAL_MODEL]]
  stopifnot(inherits(fit, "cv.glmnet"), length(scores) == 8L,
            "TCGA" %in% names(scores), all(names(scores) %in% names(inputs)))
  b <- as.matrix(coef(fit, s = fit$lambda.min))
  beta <- setNames(as.numeric(b[, 1]), rownames(b))
  beta <- beta[names(beta) != "(Intercept)"]
  selected <- names(beta)[beta != 0]
  stopifnot(setequal(selected, SIGNATURE_TFS))
  prepared <- lapply(names(scores), function(ds) {
    d <- as.data.frame(inputs[[ds]], check.names = FALSE)
    names(d) <- gsub("-", ".", names(d), fixed = TRUE)
    rs <- as.data.frame(scores[[ds]])
    stopifnot(!anyDuplicated(names(d)), all(c("ID", names(beta)) %in% names(d)),
      all(c("ID", "RS") %in% names(rs)), !anyNA(d$ID), !anyNA(rs$ID),
      !anyDuplicated(d$ID), !anyDuplicated(rs$ID))
    ix <- match(rs$ID, d$ID); stopifnot(!anyNA(ix))
    d <- d[ix, , drop = FALSE]
    X <- vapply(names(beta), function(g) {
      v <- suppressWarnings(as.numeric(as.character(d[[g]])))
      stopifnot(!any(is.infinite(v)), !all(is.na(v)))
      v[is.na(v)] <- mean(v, na.rm = TRUE); v
    }, numeric(nrow(d)))
    rownames(X) <- rs$ID
    pred <- as.numeric(predict(fit, newx = X, s = fit$lambda.min, type = "link"))
    original <- as.numeric(as.character(rs$RS))
    stopifnot(all(is.finite(original)))
    tolerance <- 1e-7 * max(1, max(abs(original)))
    error <- max(abs(pred - original))
    if (error > tolerance) stop(ds, ": saved risk scores do not match original prediction inputs")
    list(X = X, prediction = pred, tolerance = tolerance, error = error)
  })
  names(prepared) <- names(scores)
  mean_tcga <- colMeans(prepared$TCGA$X)
  baseline <- sum(mean_tcga * beta)
  tables <- file.path(output_dir, "tables"); figures <- file.path(output_dir, "figures")
  for (d in c(tables, figures)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  all_shap <- list(); importance <- list(); audit <- list(); dots <- list()
  for (ds in names(prepared)) {
    z <- prepared[[ds]]
    S <- sweep(sweep(z$X, 2, mean_tcga, "-"), 2, beta, "*")
    error <- max(abs(baseline + rowSums(S) - z$prediction))
    if (error > z$tolerance) stop(ds, ": SHAP additivity failed")
    all_shap[[ds]] <- S
    write.csv(data.frame(Sample = rownames(S), S, check.names = FALSE),
      file.path(tables, paste0(ds, "_SHAP_values.csv")), row.names = FALSE)
    importance[[ds]] <- data.frame(Dataset = ds, Gene = selected,
      Mean_abs_SHAP = colMeans(abs(S[, selected, drop = FALSE])), Coefficient = beta[selected])
    audit[[ds]] <- data.frame(Dataset = ds, N = nrow(S), Risk_error = z$error, Additivity_error = error)
    dots[[ds]] <- do.call(rbind, lapply(selected, function(g) {
      x <- z$X[, g]; span <- diff(range(x))
      data.frame(Dataset = ds, Gene = g, SHAP = S[, g], Expression = if (span == 0) .5 else (x - min(x))/span)
    }))
  }
  imp <- do.call(rbind, importance); audits <- do.call(rbind, audit)
  write.csv(imp, file.path(tables, "All_cohorts_SHAP_importance.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = names(beta), Coefficient = beta), file.path(tables, "Model_coefficients.csv"), row.names = FALSE)
  write.csv(audits, file.path(tables, "Prediction_and_additivity_audit.csv"), row.names = FALSE)
  dat <- do.call(rbind, dots)
  order <- imp$Gene[imp$Dataset == "TCGA"][order(imp$Mean_abs_SHAP[imp$Dataset == "TCGA"])]
  dat$Gene <- factor(dat$Gene, levels = order)
  p <- ggplot(dat, aes(SHAP, Gene, color = Expression)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(position = position_jitter(width = 0, height = .12, seed = 520), size = .6, alpha = .8) +
    facet_wrap(~Dataset, ncol = 4) + scale_color_gradient(low = "#2166AC", high = "#B2182B") +
    labs(x = "SHAP value (Cox log relative risk)", y = NULL) + theme_bw()
  ggsave(file.path(figures, "All_cohorts_SHAP_scatter.pdf"), p, width = 18, height = 9)
  saveRDS(list(shap = all_shap, background_mean = mean_tcga, baseline = baseline, audit = audits),
    file.path(output_dir, "exact_shap.rds"))
  invisible(list(shap = all_shap, audit = audits))
}
if (!interactive()) run_analysis()
