test_that("centrality returns a tidy one-row-per-node data frame", {
  S <- 0.4^abs(outer(1:6, 1:6, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  ct <- net_centralities(fit)
  expect_s3_class(ct, "data.frame")
  expect_named(ct, c("node", "strength", "expected_influence"))
  expect_equal(nrow(ct), 6)
  expect_equal(ct$node, fit$labels)
})

test_that("strength and expected_influence match direct row sums", {
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 300)
  g <- fit$graph
  ct <- net_centralities(fit)
  expect_equal(ct$strength, unname(rowSums(abs(g))))
  expect_equal(ct$expected_influence, unname(rowSums(g)))
})

test_that("centrality accepts a bare matrix", {
  m <- matrix(c(0, 0.3, -0.2, 0.3, 0, 0.1, -0.2, 0.1, 0), 3, 3)
  colnames(m) <- c("a", "b", "c")
  ct <- net_centralities(m)
  expect_equal(ct$node, c("a", "b", "c"))
  expect_equal(ct$strength, c(0.5, 0.4, 0.3))
})
