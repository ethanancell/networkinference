test_that("Expected results for Gaussian", {
  set.seed(1)
  A <- matrix(rnorm(n = 100^2, mean = 0, sd = 5), nrow = 100)

  # Split A
  A_split <- split_matrix(A, 'gaussian', 0.3, tau = 5)
  Atr <- A_split$Atr
  Ate <- A_split$Ate

  # "Estimate" communities
  z_hat <- rep(1:4, each = 25)
  u_mat <- matrix(c(1, -1, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, -1), nrow = 4)
  u <- as.vector(u_mat)

  gaussian_results <- infer_network(Ate, u = u, communities = z_hat,
                                    distribution = 'gaussian', epsilon = 0.3,
                                    K = 4, Atr = Atr, tau = 5)

  # Make sure that a 99.99% confidence interval contains 0 (:
  # (I guess this means this unit test will fail 0.01% of the time)
  moe <- sqrt(gaussian_results$estimate_variance) * qnorm(0.99995)
  upper_bound <- gaussian_results$estimate + moe
  lower_bound <- gaussian_results$estimate - moe

  expect_true(0 >= lower_bound)
  expect_true(0 <= upper_bound)
})
