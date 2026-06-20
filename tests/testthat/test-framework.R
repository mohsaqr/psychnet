mk <- function(seed, n = 150, p = 5) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(0.4^abs(outer(1:p, 1:p, "-")))
  colnames(X) <- paste0("V", seq_len(p))
  X
}

test_that("nct returns valid invariants and p-values", {
  fit <- nct(mk(1), mk(2), iter = 30)
  expect_s3_class(fit, "psychnet_nct")
  expect_true(fit$M$p_value >= 0 && fit$M$p_value <= 1)
  expect_true(fit$S$p_value >= 0 && fit$S$p_value <= 1)
  expect_equal(dim(fit$nw1), c(5L, 5L))
  expect_length(fit$M$perm, 30)
})

test_that("nearest-correlation projection returns a valid correlation matrix", {
  set.seed(9)
  S <- stats::cor(matrix(stats::rnorm(6 * 8), 6, 8))   # n < p -> not PD
  P <- .nearest_pd_cor(S)
  expect_true(all(abs(diag(P) - 1) < 1e-8))
  expect_gt(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("bootstrap_network returns tidy edge and centrality CIs", {
  bs <- bootstrap_network(mk(3), n_boot = 40)
  expect_s3_class(bs, "psychnet_bootstrap")
  expect_named(bs$edges,
               c("from", "to", "observed", "mean", "lower", "upper",
                 "prop_nonzero"))
  expect_true(all(bs$edges$lower <= bs$edges$upper))
  expect_true(all(bs$edges$prop_nonzero >= 0 & bs$edges$prop_nonzero <= 1))
  expect_equal(nrow(bs$edges), 10)              # 5 nodes -> 10 upper-tri edges
  expect_true(all(bs$centrality$strength_lower <= bs$centrality$strength_upper))
})

test_that("centrality_stability returns CS-coefficients in [0,1]", {
  cs <- centrality_stability(mk(4), drop_prop = c(0.3, 0.5, 0.7), iter = 15)
  expect_s3_class(cs, "psychnet_stability")
  expect_true(all(cs$cs >= 0 & cs$cs <= 1))
  expect_named(cs$table,
               c("measure", "drop_prop", "mean_cor", "prop_above"))
  # larger drop proportions never increase stability
  str_tab <- cs$table[cs$table$measure == "strength", ]
  expect_true(str_tab$mean_cor[1] + 1e-9 >= str_tab$mean_cor[nrow(str_tab)])
})
