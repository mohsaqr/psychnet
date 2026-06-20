# Nonparametric bootstrap for edge-weight and centrality accuracy, adapted from
# Nestimate's boot_glasso into a general resample-and-re-estimate loop over any
# psychnet estimator (Epskamp, Borsboom & Fried 2018).

#' Bootstrap a psychometric network
#'
#' Resamples observations with replacement, re-estimates the network on each
#' resample, and summarizes the sampling distribution of every edge weight and
#' node centrality (mean, percentile confidence interval, and edge inclusion
#' proportion).
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator (see [estimate_network()]). Default `"EBICglasso"`.
#' @param n_boot Number of bootstrap resamples. Default 1000.
#' @param ci Confidence level for percentile intervals. Default 0.95.
#' @param labels Optional node labels.
#' @param ... Passed to the estimator.
#' @return An object of class `psychnet_bootstrap` with tidy `$edges` and
#'   `$centrality` data frames and the observed network in `$observed`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- bootstrap_network(x, n_boot = 50)
#' head(bs$edges)
#' @export
bootstrap_network <- function(data, method = "EBICglasso", n_boot = 1000L,
                              ci = 0.95, labels = NULL, ...) {
  stopifnot(n_boot >= 1L, ci > 0, ci < 1)
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  if (is.null(labels)) labels <- colnames(mat)

  obs <- estimate_network(mat, method = method, labels = labels, ...)
  p <- obs$n_nodes
  ut <- upper.tri(obs$graph)
  obs_edges <- obs$graph[ut]
  obs_cent  <- centrality(obs)

  alpha <- (1 - ci) / 2
  edge_boot <- matrix(NA_real_, n_boot, length(obs_edges))
  str_boot  <- matrix(NA_real_, n_boot, p)
  ei_boot   <- matrix(NA_real_, n_boot, p)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    fit <- tryCatch(
      estimate_network(mat[idx, , drop = FALSE], method = method,
                       labels = labels, ...),
      error = function(e) NULL)
    if (is.null(fit)) next
    edge_boot[b, ] <- fit$graph[ut]
    ct <- centrality(fit)
    str_boot[b, ] <- ct$strength
    ei_boot[b, ]  <- ct$expected_influence
  }

  qci <- function(v) stats::quantile(v, c(alpha, 1 - alpha), na.rm = TRUE,
                                     names = FALSE)
  ij <- which(ut, arr.ind = TRUE)
  edge_ci <- t(apply(edge_boot, 2L, qci))
  edges <- data.frame(
    from = labels[ij[, 1L]], to = labels[ij[, 2L]],
    observed = obs_edges,
    mean = colMeans(edge_boot, na.rm = TRUE),
    lower = edge_ci[, 1L], upper = edge_ci[, 2L],
    prop_nonzero = colMeans(abs(edge_boot) > 1e-12, na.rm = TRUE),
    stringsAsFactors = FALSE, row.names = NULL)

  str_ci <- t(apply(str_boot, 2L, qci))
  ei_ci  <- t(apply(ei_boot, 2L, qci))
  cent <- data.frame(
    node = labels,
    strength = obs_cent$strength,
    strength_lower = str_ci[, 1L], strength_upper = str_ci[, 2L],
    expected_influence = obs_cent$expected_influence,
    ei_lower = ei_ci[, 1L], ei_upper = ei_ci[, 2L],
    stringsAsFactors = FALSE, row.names = NULL)

  structure(list(observed = obs, edges = edges, centrality = cent,
                 n_boot = n_boot, ci = ci, method = obs$method),
            class = "psychnet_bootstrap")
}

#' Print a network bootstrap
#'
#' @param x A `psychnet_bootstrap` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_bootstrap <- function(x, ...) {
  cat(sprintf("<psychnet_bootstrap> %s, %d resamples, %.0f%% CI\n",
              x$method, x$n_boot, 100 * x$ci))
  cat(sprintf("  %d edges, %d nodes\n", nrow(x$edges), nrow(x$centrality)))
  invisible(x)
}
