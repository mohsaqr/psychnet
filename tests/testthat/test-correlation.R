test_that("cor_network reproduces stats::cor with a zero diagonal", {
  set.seed(2)
  X <- matrix(stats::rnorm(300 * 5), 300, 5)
  colnames(X) <- paste0("V", 1:5)
  net <- cor_network(X)
  ref <- stats::cor(X); diag(ref) <- 0
  expect_equal(unname(net$graph), unname(ref))
  expect_false(net$directed)
})

test_that("pcor_network matches the inverse-correlation formula", {
  set.seed(3)
  X <- matrix(stats::rnorm(400 * 4), 400, 4)
  net <- pcor_network(X)
  S <- stats::cor(X)
  expected <- -stats::cov2cor(solve(S)); diag(expected) <- 0
  expect_equal(unname(net$graph), unname(expected), tolerance = 1e-10)
})

test_that("as.data.frame returns a tidy edge list", {
  set.seed(4)
  X <- matrix(stats::rnorm(200 * 4), 200, 4)
  net <- cor_network(X, threshold = 0.05)
  ed <- as.data.frame(net)
  expect_named(ed, c("from", "to", "weight"))
  expect_true(all(abs(ed$weight) >= 0.05))
  expect_equal(nrow(ed), net$n_edges)
})

test_that("estimate_network dispatches and aliases resolve", {
  set.seed(5)
  X <- matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
  expect_equal(estimate_network(X, "cor")$method, "cor")
  expect_equal(estimate_network(X, "pcor")$method, "pcor")
  expect_equal(estimate_network(X, "glasso")$method, "EBICglasso")
  expect_equal(estimate_network(X, "EBICglasso")$method, "EBICglasso")
})
