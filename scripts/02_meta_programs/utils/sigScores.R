# ============================================================================
# Script: sigScores.R
# Module: 02_meta_programs/utils
# Description:
#   Signature scoring utilities including gene set filtering, base scores,
#   and expression-bin-corrected signature scores.
# ============================================================================

#' Filter genes in signatures according to reference
#'
#' @param sigs A character vector of genes or list of character vectors to filter.
#' @param ref Reference genes to filter sigs according to.
#' @param conserved The minimum allowed fraction of genes retained after
#'   filtering. Default 0.7.
#' @return A filtered list of sigs. Returns NULL if no genes remain.
#'   sigs with less than conserved fraction of genes retained will be removed.
#'   Returns input sigs if no genes are missing.
#' @rdname filter_sigs
#' @export
filter_sigs <- function(sigs, ref, conserved = 0.7) {
    if (all(unlist(sigs) %in% ref)) {
        return(sigs)
    }

    ngr0 <- length(sigs)
    nge0 <- lengths(sigs)

    sigs <- sapply(sigs, function(group) group[group %in% ref], simplify = FALSE)
    nge1 <- lengths(sigs)

    if (all(nge1 == 0)) {
        return(NULL)
    }

    frac.conserved <- nge1 / nge0
    sigs <- sigs[frac.conserved >= conserved]

    frac.conserved <- round(frac.conserved, 2)
    frac.conserved <- paste(names(frac.conserved), frac.conserved, sep = ": ")

    ngr1 <- length(sigs)
    if (ngr1 == 0) {
        stop("No sigs left to score with.")
    }

    nge1 <- nge1[names(sigs)]
    nge0 <- nge0[names(sigs)]

    if (ngr1 < ngr0) {
        warning("Removed ", ngr0 - ngr1, " out of ", ngr0,
                " sigs with < ", conserved * 100,
                "% genes retained after filtering...")
    }

    if (any(nge1 != nge0)) {
        warning("Some genes were filtered out. Fractions retained:", "\n",
                paste0(frac.conserved, collapse = "\n"))
    }

    sigs
}

#' Basic Scoring of Matrix by Gene Signatures
#'
#' Average expression level of each column in m for each signature in sigs.
#'
#' @param m A non-centered matrix of genes X cells/samples.
#' @param sigs A character vector of genes or list of character vectors.
#' @param conserved.genes Minimum fraction of genes retained after filtering
#'   that is allowed. Default 0.7.
#' @return A data frame of scores.
#' @rdname baseScores
#' @export
baseScores <- function(m, sigs, conserved.genes = 0.7) {
    if (is.character(sigs)) sigs <- list(sigs)
    sigs <- filter_sigs(sigs, ref = rownames(m), conserved = conserved.genes)
    sapply(sigs, function(sig) colMeans(m[sig, , drop = FALSE]))
}


#' Score a Matrix by Gene Signatures
#'
#' Score a matrix by gene signatures with options for expression-bin-matched
#' normalization. Raw scores are generated with baseScores (average expression
#' of genes in the signature). Options to normalize by expression-bin-matched
#' controls or by centering.
#'
#' @param m A non-centered matrix of genes X cells/samples.
#' @param sigs A character vector of genes or list of character vectors.
#' @param groups If scores should be calculated intra-tumour, a list of cell
#'   IDs by sample. Default NULL.
#' @param center.rows Logical; row-center m before scoring. Default TRUE.
#' @param center Logical; center scores by mean expression. Default TRUE.
#' @param expr.center Logical; normalize by expression-bin-matched sigs.
#'   Takes precedence over center if both TRUE. Default TRUE.
#' @param expr.bin.m Matrix used to bin genes. If NULL, uses m. Default NULL.
#' @param expr.bins Pre-computed bin vector. Default NULL.
#' @param expr.sigs Pre-computed expression control sigs. Default NULL.
#' @param expr.nbin Number of expression bins. Default 30.
#' @param expr.binsize Number of bin-matched genes per gene. Default 100.
#' @param conserved.genes Minimum fraction of genes retained. Default 0.7.
#' @param replace Logical; allow replacement in bin sampling. Default FALSE.
#' @return Data frame of cell/sample scores.
#' @rdname sigScores
#' @export
sigScores <- function(m,
                     sigs,
                     groups = NULL,
                     center.rows = TRUE,
                     center = TRUE,
                     expr.center = TRUE,
                     expr.bin.m = NULL,
                     expr.bins = NULL,
                     expr.sigs = NULL,
                     expr.nbin = 30,
                     expr.binsize = 100,
                     conserved.genes = 0.7,
                     replace = FALSE) {

    if (is.character(sigs)) sigs <- list(sigs)
    sigs <- filter_sigs(sigs, ref = rownames(m), conserved = conserved.genes)

    Args <- mget(ls())

    if (!is.null(groups)) {
        Args$expr.bin.m <- m
        Args <- Args[names(Args) != "m"]
        mlist <- sapply(groups, function(sample) m[, sample], simplify = FALSE)
        res <- sapply(mlist, function(m) {
            do.call(.sigScores, c(list(m = m), Args))
        }, simplify = FALSE)

        cells <- unlist(sapply(res, rownames, simplify = FALSE))
        res <- do.call(rbind.data.frame, res)
        rownames(res) <- cells
    } else {
        res <- do.call(.sigScores, Args)
    }

    res
}


#' Internal signature scoring function
#'
#' @inheritParams sigScores
#' @keywords internal
.sigScores <- function(m,
                      sigs,
                      groups = NULL,
                      center.rows = TRUE,
                      center = TRUE,
                      expr.center = TRUE,
                      expr.bin.m = NULL,
                      expr.bins = NULL,
                      expr.sigs = NULL,
                      expr.nbin = 30,
                      expr.binsize = 100,
                      conserved.genes = 0.7,
                      replace = FALSE) {

    if (center.rows) {
        scores <- baseScores(m = rowcenter(m), sigs = sigs, conserved.genes = conserved.genes)
    } else {
        scores <- baseScores(m = m, sigs = sigs, conserved.genes = conserved.genes)
    }
    if (!center) expr.center <- FALSE

    if (expr.center) {
        if (is.null(expr.sigs)) {
            if (is.null(expr.bins)) {
                if (is.null(expr.bin.m)) expr.bin.m <- m
                expr.bins <- bin(expr.bin.m, breaks = expr.nbin)
                stopifnot(all(unlist(sigs) %in% names(expr.bins)))
            }
            expr.sigs <- sapply(sigs,
                               binmatch,
                               bins = expr.bins,
                               n = expr.binsize,
                               replace = replace,
                               simplify = FALSE)
            names(expr.sigs) <- names(sigs)
        }

        if (center.rows) {
            expr.scores <- baseScores(m = rowcenter(m), sigs = expr.sigs)
        } else {
            expr.scores <- baseScores(m = m, sigs = expr.sigs)
        }
        scores <- scores - expr.scores
    } else if (center) {
        if (!is.null(expr.bin.m)) center.scores <- colMeans(expr.bin.m)
        else center.scores <- colMeans(m)
        scores <- sweep(scores, MARGIN = 1, STATS = center.scores, FUN = "-")
    }

    rows <- rownames(scores)
    scores <- as.data.frame(scores)
    rownames(scores) <- rows
    scores
}

#' Score a Matrix with Marker Gene Sets Of Normal Cell Types
#'
#' Wrapper around sigScores using normal cell type markers.
#'
#' @param m A non-centered matrix of genes X cells/samples.
#' @param ... Other arguments passed to sigScores.
#' @return Data frame of cell/sample scores.
#' @rdname markerScores
#' @export
markerScores <- function(m, ...) {
    sigs <- suppressWarnings(filter_sigs(Markers_Normal, ref = rownames(m), conserved = 0.4))
    sigs <- sigs[lengths(sigs) >= 2]
    sigScores(m = m, sigs = sigs, conserved.genes = conserved, ...)
}
