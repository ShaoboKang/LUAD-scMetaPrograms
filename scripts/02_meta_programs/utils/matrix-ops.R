# ============================================================================
# Script: matrix-ops.R
# Module: 02_meta_programs/utils
# Description:
#   Matrix operation utilities including centering, detection counting,
#   and dimension checks.
# ============================================================================

#' Center a matrix column-wise
#'
#' @param m A matrix or Matrix.
#' @param by Either "mean", "median" or a numeric vector of length equal to
#'   the number of columns of m. Default "mean".
#' @return Column-centered matrix.
#' @rdname colcenter
#' @export
colcenter <- function(m, by = "mean") {
    m <- as.matrix(m)
    if (by == "mean") by <- colMeans(m, na.rm = TRUE)
    else if (by == "median") by <- matrixStats::colMedians(m, na.rm = TRUE)
    else stopifnot(is.numeric(by) & length(by) == ncol(m))
    scale(m, center = by, scale = FALSE)
}

#' Center a matrix row-wise
#'
#' @param m A matrix or Matrix.
#' @param by Either "mean", "median" or a numeric vector of length equal to
#'   the number of rows of m. Default "mean".
#' @return Row-centered matrix.
#' @rdname rowcenter
#' @export
rowcenter <- function(m, by = "mean") {
    m <- as.matrix(m)
    if (by == "mean") by <- rowMeans(m, na.rm = TRUE)
    else if (by == "median") by <- matrixStats::rowMedians(m, na.rm = TRUE)
    else stopifnot(is.numeric(by) & length(by) == nrow(m))
    t(scale(t(m), center = by, scale = FALSE))
}

#' Number of non-zero values per column
#'
#' @param m A matrix.
#' @param value Value to compare against. Default 0.
#' @param method Comparison method. Default "notequal".
#' @param counts If TRUE, return counts; if FALSE, return logical matrix.
#' @return Numeric vector of counts per column.
#' @rdname coldetected
#' @export
coldetected <- function(m, value = 0,
                       method = c("notequal", "greaterthan", "lessthan", "equal"),
                       counts = TRUE) {
    method <- match.arg(method)
    m <- as.matrix(m)
    if (method == "notequal") {
        bool <- m != value
        if (!counts) return(bool)
        res <- matrixStats::colCounts(bool)
    }
    if (method == "equal") {
        bool <- m == value
        if (!counts) return(bool)
        res <- matrixStats::colCounts(bool)
    }
    if (method == "greaterthan") {
        bool <- m > value
        if (!counts) return(bool)
        res <- matrixStats::colCounts(bool)
    }
    if (method == "lessthan") {
        bool <- m < value
        if (!counts) return(bool)
        res <- matrixStats::colCounts(bool)
    }
    stats::setNames(res, colnames(m))
}

#' Number of non-zero values per row
#'
#' @param m A matrix.
#' @param value Value to compare against. Default 0.
#' @param method Comparison method. Default "notequal".
#' @param counts If TRUE, return counts; if FALSE, return logical matrix.
#' @return Numeric vector of counts per row.
#' @rdname rowdetected
#' @export
rowdetected <- function(m, value = 0,
                       method = c("notequal", "greaterthan", "lessthan", "equal"),
                       counts = TRUE) {
    method <- match.arg(method)
    m <- as.matrix(m)
    if (method == "notequal") {
        bool <- m != value
        if (!counts) return(bool)
        res <- matrixStats::rowCounts(bool)
    }
    if (method == "equal") {
        bool <- m == value
        if (!counts) return(bool)
        res <- matrixStats::rowCounts(bool)
    }
    if (method == "greaterthan") {
        bool <- m > value
        if (!counts) return(bool)
        res <- matrixStats::rowCounts(bool)
    }
    if (method == "lessthan") {
        bool <- m < value
        if (!counts) return(bool)
        res <- matrixStats::rowCounts(bool)
    }
    stats::setNames(res, rownames(m))
}

#' Dimensions for many matrices
#'
#' @param mats A list of matrices (or a single matrix).
#' @return dim for each matrix provided.
#' @rdname dims
#' @export
dims <- function(mats) {
    if (!is.null(dim(mats))) {
        return(dim(mats))
    }
    sapply(mats, dim, simplify = TRUE)
}

#' Number of columns for many matrices
#'
#' @param mats A list of matrices (or a single matrix).
#' @return ncol for each matrix provided.
#' @rdname ncols
#' @export
ncols <- function(mats) {
    if (!is.null(dim(mats))) {
        return(ncol(mats))
    }
    sapply(mats, ncol, simplify = TRUE)
}

#' Number of rows for many matrices
#'
#' @param mats A list of matrices (or a single matrix).
#' @return nrow for each matrix provided.
#' @rdname nrows
#' @export
nrows <- function(mats) {
    if (!is.null(dim(mats))) {
        return(nrow(mats))
    }
    sapply(mats, nrow, simplify = TRUE)
}

#' Check if object has dimensions
#'
#' @param x Any R object.
#' @return Logical.
#' @export
has_dim <- function(x) {
    if (is.data.frame(x)) x <- as.matrix(x)
    !is.null(attr(x, "dim"))
}

#' Split a matrix by column names
#'
#' @param m A matrix.
#' @param by Character vector of column names to extract.
#' @return List with x (selected columns) and y (remaining columns).
#' @export
split_matrix <- function(m, by) {
    stopifnot(has_dim(m))
    stopifnot(is.character(by))
    stopifnot(all(by %in% colnames(m)))
    list(x = m[, by, drop = FALSE], y = m[, !colnames(m) %in% by, drop = FALSE])
}

#' Check equal number of rows
#'
#' @param m1 First matrix.
#' @param m2 Second matrix.
#' @return Logical.
#' @export
have_equal_nrows <- function(m1, m2) {
    nrow(m1) == nrow(m2)
}

#' Check equal row names
#'
#' @param m1 First matrix.
#' @param m2 Second matrix.
#' @return Logical.
#' @export
have_equal_rownames <- function(m1, m2) {
    all(rownames(m1) == rownames(m2))
}

#' Check if matrix is square
#'
#' @param m A matrix.
#' @return Logical.
#' @export
is_square <- function(m) {
    nrow(m) == ncol(m)
}

#' Check equal dimensions
#'
#' @param m1 First matrix.
#' @param m2 Second matrix.
#' @return Logical.
#' @export
have_equal_dims <- function(m1, m2) {
    identical(dim(m1), dim(m2))
}

#' Check if matrix is a correlation matrix
#'
#' @param m A matrix.
#' @return Logical.
#' @export
is_cor <- function(m) {
    rg <- range(m)
    if ((is_square(m)) & (rg[1] >= -1) & (rg[2] <= 1)) {
        dg <- unique(diag(m))
        return(length(dg) == 1 & dg == 1)
    }
    FALSE
}

#' Check if matrix is symmetric
#'
#' @param m A matrix.
#' @return Logical.
#' @export
is_symm <- function(m) {
    (is_square(m)) && (sum(m == t(m)) == nrow(m)^2)
}
