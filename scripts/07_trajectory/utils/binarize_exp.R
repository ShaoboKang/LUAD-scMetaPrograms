# ============================================================================
# Utility: binarize_exp.R
# Description: Binarize single-cell expression data using mixture modeling
#   or a fixed cutoff.
# ============================================================================

#' Binarize single-cell expression
#'
#' Fits a two-component Gaussian mixture per gene to determine an on/off
#'.cutoff, or uses a fixed cutoff if requested.
#'
#' @param sce A SingleCellExperiment with an 'expdata' assay.
#' @param fix_cutoff Logical; if TRUE use a fixed cutoff.
#' @param binarize_cutoff Numeric cutoff when fix_cutoff is TRUE.
#' @param ncores Number of cores for parallel processing.
#' @return Modified SingleCellExperiment with binary assay and rowData.
#' @importFrom parallel mclapply
binarize_exp <- function(sce, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 3) {
  zerop_g <- c()
  expdata <- assays(sce)$expdata
  for (i in seq_len(nrow(expdata))) {
    zp <- length(which(expdata[i, ] == 0)) / ncol(expdata)
    zerop_g <- c(zerop_g, zp)
  }

  if (fix_cutoff == TRUE) {
    expdata <- assays(sce)$expdata
    is.na(expdata) <- assays(sce)$expdata == 0
    exp_reduced_binary <- as.matrix((expdata > binarize_cutoff) + 0)
    exp_reduced_binary[is.na(exp_reduced_binary)] <- 0
    assays(sce)$binary <- exp_reduced_binary
    oup_binary <- data.frame(
      geneID = rownames(sce),
      zerop_gene = zerop_g,
      passBinary = TRUE
    )
    rowData(sce) <- oup_binary
  } else {
    expdata <- assays(sce)$expdata
    log_counts_add <- expdata + matrix(
      rnorm(nrow(expdata) * ncol(expdata), mean = 0, sd = 0.1),
      nrow(expdata), ncol(expdata)
    )
    oup_binary <- do.call(rbind, parallel::mclapply(rownames(log_counts_add), function(i_gene) {
      set.seed(42)
      tmp_mix <- mixtools::normalmixEM(log_counts_add[i_gene, ], k = 2)
      if (tmp_mix$mu[1] < tmp_mix$mu[2]) {
        tmp_oup <- data.frame(
          geneID = i_gene,
          mu1 = tmp_mix$mu[1],
          mu2 = tmp_mix$mu[2],
          sigma1 = tmp_mix$sigma[1],
          sigma2 = tmp_mix$sigma[2],
          lambda1 = tmp_mix$lambda[1],
          lambda2 = tmp_mix$lambda[2],
          loglik = tmp_mix$loglik
        )
      } else {
        tmp_oup <- data.frame(
          geneID = i_gene,
          mu1 = tmp_mix$mu[2],
          mu2 = tmp_mix$mu[1],
          sigma1 = tmp_mix$sigma[2],
          sigma2 = tmp_mix$sigma[1],
          lambda1 = tmp_mix$lambda[2],
          lambda2 = tmp_mix$lambda[1],
          loglik = tmp_mix$loglik
        )
      }
      return(tmp_oup)
    }, mc.cores = ncores))

    oup_binary$passBinary <- TRUE
    oup_binary[oup_binary$lambda1 < 0.1, ]$passBinary <- FALSE
    oup_binary[oup_binary$lambda2 < 0.1, ]$passBinary <- FALSE
    oup_binary[(oup_binary$mu2 - oup_binary$mu1) <
                 (oup_binary$sigma1 + oup_binary$sigma2), ]$passBinary <- FALSE
    oup_binary$root <- -1

    for (i_gene in oup_binary[oup_binary$passBinary == TRUE, ]$geneID) {
      tmp_mix <- oup_binary[oup_binary$geneID == i_gene, ]
      tmp_int <- stats::uniroot(function(x, l1, l2, mu1, mu2, sd1, sd2) {
        dnorm(x, m = mu1, sd = sd1) * l1 - dnorm(x, m = mu2, sd = sd2) * l2
      }, interval = c(tmp_mix$mu1, tmp_mix$mu2),
      l1 = tmp_mix$lambda1, mu1 = tmp_mix$mu1, sd1 = tmp_mix$sigma1,
      l2 = tmp_mix$lambda2, mu2 = tmp_mix$mu2, sd2 = tmp_mix$sigma2)
      oup_binary[oup_binary$geneID == i_gene, ]$root <- tmp_int$root
    }

    bin_log_counts <- as.matrix(expdata[oup_binary$geneID, ])
    bin_log_counts <- t(scale(t(bin_log_counts), scale = FALSE, center = oup_binary$root))
    bin_log_counts[bin_log_counts >= 0] <- 1
    bin_log_counts[bin_log_counts < 0] <- 0
    assays(sce)$binary <- bin_log_counts
    oup_binary$zerop_gene <- zerop_g
    rowData(sce) <- oup_binary
  }
  return(sce)
}
