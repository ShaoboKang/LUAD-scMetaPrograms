# ============================================================================
# Script: 00_tools_robust_nmf_programs.R
# Module: 02_meta_programs/utils
# Description:
#   Utility functions for selecting robust NMF programs and custom color palettes.
# ============================================================================

#' Select robust NMF programs
#'
#' Select robust nonnegative matrix factorization (NMF) programs based on
#' intra-sample and inter-sample overlap thresholds.
#'
#' @param nmf_programs A list; each element contains a matrix with NMF programs
#'   (top 50 genes) generated for a specific sample using different NMF ranks.
#' @param intra_min Minimum overlap with a program from the same sample
#'   (for selecting robust programs). Default 35.
#' @param intra_max Maximum overlap with a program from the same sample
#'   (for removing redundant programs). Default 10.
#' @param inter_filter Logical; whether programs should be filtered based on
#'   their similarity to programs of other samples. Default TRUE.
#' @param inter_min Minimum overlap with a program from another sample. Default 10.
#' @return Character vector with the selected NMF program names.
#' @export
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10,
                                inter_filter = TRUE, inter_min = 10) {
  # Select NMF programs based on minimum overlap with other programs from the same sample
  intra_intersect <- lapply(nmf_programs, function(z) {
    apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y))))
  })
  intra_intersect_max <- lapply(intra_intersect, function(x) {
    apply(x, 2, function(y) sort(y, decreasing = TRUE)[2])
  })
  nmf_sel <- lapply(names(nmf_programs), function(x) {
    nmf_programs[[x]][, intra_intersect_max[[x]] >= intra_min]
  })
  names(nmf_sel) <- names(nmf_programs)

  # Select programs based on max intra-sample overlap and min inter-sample overlap
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  require(parallel)
  cl <- makeCluster(8, type = "FORK")
  inter_intersect <- parApply(cl = cl, nmf_sel_unlist, 2,
                              function(x) apply(nmf_sel_unlist, 2, function(y) length(intersect(x, y))))
  stopCluster(cl)

  final_filter <- NULL
  for (i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = TRUE),
                         grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = TRUE)
    if (inter_filter == TRUE) {
      b <- b[b >= inter_min]
    }
    if (length(b) > 1) {
      c <- names(b[1])
      for (y in 2:length(b)) {
        if (max(inter_intersect[c, names(b[y])]) <= intra_max) {
          c <- c(c, names(b[y]))
        }
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)
}

# ---- Custom color palette ----

#' Custom magma color palette
#'
#' A extended magma-like color palette for heatmaps.
#'
#' @export
custom_magma <- c(
  colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10),
  rev(viridis::magma(323, begin = 0.18))
)
