# Correctness is certified against the convex objective itself (KKT residual),
# so these tests need no reference solver.

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
compound <- function(p, rho) { m <- matrix(rho, p, p); diag(m) <- 1; m }

test_that("ebic_glasso returns the certified global optimum", {
  for (S in list(ar1(6, 0.5), ar1(8, 0.6), compound(5, 0.4))) {
    fit <- ebic_glasso(cor_matrix = S, n = 250)
    expect_s3_class(fit, "psychnet")
    expect_lt(fit$kkt, 1e-7)
    expect_equal(glasso_kkt(fit$precision, S, fit$lambda), fit$kkt)
  }
})

test_that("glasso_kkt flags a non-optimal precision matrix", {
  S <- ar1(5, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 250)
  bad <- fit$precision
  bad[1, 2] <- bad[2, 1] <- bad[1, 2] + 0.3
  expect_gt(glasso_kkt(bad, S, fit$lambda),
            glasso_kkt(fit$precision, S, fit$lambda))
})

test_that("ebic_glasso graph is a symmetric partial-correlation matrix", {
  S <- ar1(7, 0.5)
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  g <- fit$graph
  expect_equal(g, t(g))
  expect_true(all(diag(g) == 0))
  expect_true(all(abs(g) <= 1 + 1e-8))
})

test_that("ebic_glasso runs from raw data and respects threshold", {
  set.seed(1)
  p <- 6
  X <- matrix(stats::rnorm(400 * p), 400, p) %*% chol(ar1(p, 0.5))
  colnames(X) <- paste0("V", seq_len(p))
  fit  <- ebic_glasso(as.data.frame(X), gamma = 0.5)
  fitT <- ebic_glasso(as.data.frame(X), gamma = 0.5, threshold = 0.05)
  expect_lt(fit$kkt, 1e-6)
  expect_true(all(abs(fitT$graph[fitT$graph != 0]) >= 0.05))
  expect_lte(fitT$n_edges, fit$n_edges)
})
