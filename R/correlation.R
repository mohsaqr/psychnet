# Unregularized correlation and partial-correlation networks (pure base R).

# Coerce a data frame / matrix to a clean numeric matrix: keep numeric columns
# with non-zero variance, drop incomplete rows, guarantee column names.
#' @noRd
.as_numeric_matrix <- function(data) {
  if (is.null(data)) stop("`data` is required.", call. = FALSE)
  if (is.matrix(data)) {
    mat <- data
    storage.mode(mat) <- "double"
  } else if (is.data.frame(data)) {
    num <- vapply(data, is.numeric, logical(1))
    if (!any(num)) stop("No numeric columns in `data`.", call. = FALSE)
    mat <- as.matrix(data[, num, drop = FALSE])
  } else {
    stop("`data` must be a data frame or numeric matrix.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  vars <- apply(mat, 2L, stats::sd)
  if (any(vars == 0 | is.na(vars))) {
    mat <- mat[, vars > 0 & !is.na(vars), drop = FALSE]
  }
  if (ncol(mat) < 2L) stop("Need at least 2 usable variables.", call. = FALSE)
  mat
}

#' Correlation network
#'
#' Marginal (zero-order) association network: the Pearson correlation matrix
#' with the diagonal removed. Equivalent to `bootnet`'s `"cor"` default.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Correlation method: `"pearson"` (default), `"spearman"`, or
#'   `"kendall"`.
#' @param threshold Correlations with absolute value below this are set to zero.
#'   Default 0.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the thresholded correlation
#'   matrix.
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' cor_network(x)
#' @export
cor_network <- function(data, method = c("pearson", "spearman", "kendall"),
                        threshold = 0, labels = NULL) {
  method <- match.arg(method)
  mat <- .as_numeric_matrix(data)
  if (is.null(labels)) labels <- colnames(mat)
  g <- stats::cor(mat, method = method)
  diag(g) <- 0
  g[abs(g) < threshold] <- 0
  .new_psychnet(g, labels, method = "cor", directed = FALSE,
                n_obs = nrow(mat),
                extra = list(cor_matrix = stats::cor(mat, method = method)))
}

#' Partial correlation network
#'
#' Conditional (full-order) association network: each edge is the correlation
#' between two variables with all others partialled out, obtained from the
#' inverse correlation matrix. Equivalent to `bootnet`'s `"pcor"` default.
#'
#' @inheritParams cor_network
#' @return A `psychnet` object whose `$graph` is the thresholded
#'   partial-correlation matrix, with `$precision` and `$cor_matrix`.
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' pcor_network(x)
#' @export
pcor_network <- function(data, method = c("pearson", "spearman", "kendall"),
                         threshold = 0, labels = NULL) {
  method <- match.arg(method)
  mat <- .as_numeric_matrix(data)
  if (is.null(labels)) labels <- colnames(mat)
  S  <- stats::cor(mat, method = method)
  wi <- solve(S)
  g  <- .precision_to_pcor(wi)
  g[abs(g) < threshold] <- 0
  dimnames(g) <- list(labels, labels)
  .new_psychnet(g, labels, method = "pcor", directed = FALSE,
                n_obs = nrow(mat),
                extra = list(precision = wi, cor_matrix = S))
}
