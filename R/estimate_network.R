# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...).

#' Estimate a psychometric network
#'
#' Single entry point that routes to the requested estimator and returns a
#' common [psychnet] object, so callers can swap estimators without rewiring
#' downstream code. Mirrors `bootnet::estimateNetwork(data, default = ...)`.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method One of `"cor"`, `"pcor"`, `"EBICglasso"`, `"ising"`, `"mgm"`.
#'   Aliases: `"glasso"` -> `"EBICglasso"`, `"IsingFit"` -> `"ising"`.
#' @param threshold Absolute-weight threshold below which edges are zeroed.
#' @param gamma EBIC hyperparameter for the regularized methods.
#' @param labels Optional node labels.
#' @param ... Passed to the underlying estimator.
#' @return A `psychnet` object.
#' @examples
#' x <- matrix(stats::rnorm(200 * 5), 200, 5)
#' estimate_network(x, method = "EBICglasso")
#' estimate_network(x, method = "pcor")
#' @export
estimate_network <- function(data,
                             method = c("EBICglasso", "cor", "pcor",
                                        "ising", "mgm"),
                             threshold = 0, gamma = 0.5, labels = NULL, ...) {
  method <- .resolve_method(method)
  switch(
    method,
    cor        = cor_network(data, threshold = threshold, labels = labels, ...),
    pcor       = pcor_network(data, threshold = threshold, labels = labels, ...),
    EBICglasso = ebic_glasso(data, gamma = gamma, threshold = threshold,
                             labels = labels, ...),
    ising      = ising_fit(data, gamma = gamma, labels = labels, ...),
    mgm        = mgm_fit(data, gamma = gamma, labels = labels, ...),
    stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
}

# Resolve method name + aliases.
#' @noRd
.resolve_method <- function(method) {
  if (length(method) > 1L) method <- method[1L]
  aliases <- c(glasso = "EBICglasso", ebicglasso = "EBICglasso",
               EBICglasso = "EBICglasso", isingfit = "ising", IsingFit = "ising",
               ising = "ising", cor = "cor", correlation = "cor",
               pcor = "pcor", partial = "pcor", mgm = "mgm")
  key <- aliases[method]
  if (is.na(key)) key <- aliases[tolower(method)]
  if (is.na(key)) stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  unname(key)
}
