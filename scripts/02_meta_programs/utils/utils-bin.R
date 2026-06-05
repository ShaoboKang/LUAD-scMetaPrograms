# ============================================================================
# Script: utils-bin.R
# Module: 02_meta_programs/utils
# Description:
#   Internal validation helpers for binning and signature functions.
# ============================================================================

#' Check that argument is a single character vector
#'
#' @param arg Argument to check.
#' @keywords internal
.check_arg <- function(arg) {
    if (is.list(arg) & all(lengths(arg) >= 1)) {
        stop("Please provide a single character vector.")
    }
}

#' Check if argument contains multiple character vectors
#'
#' @param arg Argument to check.
#' @return Logical.
#' @keywords internal
.arg_is_args <- function(arg) {
    arg <- ifelse(is.list(arg) & all(lengths(arg) >= 1), TRUE, FALSE)
    return(arg)
}

#' Check sample size constraints for bin sampling
#'
#' @param bins Named list of bins.
#' @param n Requested sample size.
#' @param replace Logical; sampling with replacement.
#' @keywords internal
.check_sample_size <- function(bins, n, replace) {
    err1 <- "Error: Group is smaller than sample size <n> and <replace> == FALSE."
    err2 <- "Set <replace> = TRUE to proceed (or make <n> smaller)."
    errors <- c(err1, err2)
    if (any(lengths(bins) < n) & replace == FALSE) {
        stop(paste(errors, collapse = "\n"))
    }
}

#' Name bin IDs using values from x
#'
#' @param binIDs Numeric vector of bin IDs.
#' @param x Named numeric vector.
#' @return Named vector of bin IDs.
#' @keywords internal
.name_binIDs <- function(binIDs, x) {
    if (!is.null(names(x))) {
        names(binIDs) <- names(x)
    } else {
        names(binIDs) <- x
    }
    return(binIDs)
}

#' Check for missing values in Group relative to reference
#'
#' @param Group Character vector to check.
#' @param ref Reference vector.
#' @keywords internal
.check_missing <- function(Group, ref) {
    if (is.numeric(ref)) ref <- names(ref)
    are_missing <- !Group %in% ref
    err <- "Error: Returning missing name(s)..."
    if (any(are_missing)) {
        missing <- Group[are_missing]
        stop(paste(c(err, missing), collapse = "\n"))
    }
}

#' Check that at least one required argument is provided
#'
#' @param x Named numeric vector or NULL.
#' @param data Matrix or NULL.
#' @param bins Bin vector or NULL.
#' @keywords internal
.check_args_exist <- function(x, data, bins) {
    args <- list(x, data, bins)
    err1 <- "Not enough arguments. Provide one of:"
    err2 <- "\t1. <data> : data to be transformed to <x>"
    err3 <- "\t2. <x> : vector to be binned"
    err4 <- "\t3. <bins>: bin vector for control <Group>"
    errors <- c(err1, err2, err3, err4)
    if (all(sapply(args, is.null))) {
        stop(paste(errors, collapse = "\n"))
    }
}
