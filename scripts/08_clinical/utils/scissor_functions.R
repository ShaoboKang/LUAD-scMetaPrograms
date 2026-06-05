# ============================================================================
# Utility: scissor_functions.R
# Description: Scissor cell selection function for associating single-cell
#   states with bulk phenotypes.
# ============================================================================

#' Scissor: Identify phenotype-associated cells
#'
#' Uses network-regularized regression to select single cells that are
#' informative for a given bulk phenotype (gaussian, binomial, or cox).
#'
#' @param bulk_dataset Bulk expression matrix (genes x samples).
#' @param sc_dataset Single-cell expression matrix or Seurat object.
#' @param phenotype Phenotype vector or survival matrix.
#' @param tag Labels for phenotype classes (for gaussian/binomial).
#' @param alpha Vector of regularization parameters to search.
#' @param cutoff Maximum fraction of selected cells allowed.
#' @param family Regression family: "gaussian", "binomial", or "cox".
#' @param save_file Path to save intermediate RData.
#' @param load_file Path to load intermediate RData (skip preprocessing).
#' @return List with selected cell IDs, coefficients, and tuning parameters.
#' @importFrom Seurat GetAssayData CreateSeuratObject FindVariableFeatures
#'   ScaleData RunPCA FindNeighbors
#' @importFrom Matrix Matrix
#' @importFrom preprocessCore normalize.quantiles
Scissor <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL,
                    alpha = NULL, cutoff = 0.2,
                    family = c("gaussian", "binomial", "cox"),
                    save_file = "Scissor_inputs.RData",
                    load_file = NULL) {
  library(Seurat)
  library(Matrix)
  library(preprocessCore)

  family <- match.arg(family)

  if (is.null(load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("No common genes between the given single-cell and bulk samples.")
    }

    if (inherits(sc_dataset, "Seurat")) {
      sc_exprs <- as.matrix(GetAssayData(sc_dataset, layer = "data"))
      network <- as.matrix(sc_dataset@graphs$RNA_snn)
    } else {
      sc_exprs <- as.matrix(sc_dataset)
      seurat_tmp <- CreateSeuratObject(sc_dataset)
      seurat_tmp <- FindVariableFeatures(seurat_tmp, selection.method = "vst", verbose = FALSE)
      seurat_tmp <- ScaleData(seurat_tmp, verbose = FALSE)
      seurat_tmp <- RunPCA(seurat_tmp, features = VariableFeatures(seurat_tmp), verbose = FALSE)
      seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:10, verbose = FALSE)
      network <- as.matrix(seurat_tmp@graphs$RNA_snn)
    }

    diag(network) <- 0
    network[which(network != 0)] <- 1

    dataset0 <- as.matrix(cbind(bulk_dataset[common, ], sc_exprs[common, ]))
    dataset1 <- normalize.quantiles(dataset0)
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)

    expression_bulk <- dataset1[, seq_len(ncol(bulk_dataset))]
    expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- cor(expression_bulk, expression_cell)

    quality_check <- quantile(X)
    message("Performing quality-check for the correlations")
    message("The five-number summary of correlations: ")
    print(quality_check)
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }

    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes.")
      } else {
        message(sprintf("Current phenotype contains %d %s and %d %s samples.",
                        z[1], tag[1], z[2], tag[2]))
        message("Perform logistic regression on the given phenotypes:")
      }
    }

    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes.")
      } else {
        tmp <- paste(z, tag)
        message(paste0("Current phenotype contains ",
                       paste(tmp[seq_len(length(z) - 1)], collapse = ", "),
                       ", and ", tmp[length(z)], " samples."))
        message("Perform linear regression on the given phenotypes:")
      }
    }

    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check inputs.")
      } else {
        message("Perform cox regression on the given clinical outcomes:")
      }
    }

    save(X, Y, network, expression_bulk, expression_cell, file = save_file)
  } else {
    load(load_file)
  }

  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
               0.6, 0.7, 0.8, 0.9)
  }

  for (i in seq_along(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net",
                  alpha = alpha[i], Omega = network, nlambda = 100,
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net",
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)

    if (family == "binomial") {
      coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    } else {
      coefs <- as.numeric(fit1$Beta)
    }

    cell_pos <- colnames(X)[which(coefs > 0)]
    cell_neg <- colnames(X)[which(coefs < 0)]
    percentage <- (length(cell_pos) + length(cell_neg)) / ncol(X)

    message(sprintf("alpha = %s", alpha[i]))
    message(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.",
                    length(cell_pos), length(cell_neg)))
    message(sprintf("The percentage of selected cell is: %s%%",
                    formatC(percentage * 100, format = "f", digits = 3)))

    if (percentage < cutoff) {
      break
    }
  }

  message("Scissor selection complete.")
  return(list(
    para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family),
    Coefs = coefs,
    Scissor_pos = cell_pos,
    Scissor_neg = cell_neg
  ))
}
