# Checking that error handling is done correctly

# Gaussian tests
test_that('Gaussian splitting works', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  splitting_results1 <- split_matrix(A, 'gaussian', 0.3, sigma = 5)
  splitting_results2 <- split_matrix(A, 'gaussian', 0.5, sigma = 5)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})

# Poisson tests
test_that('Poisson splitting works', {
  A <- matrix(rpois(n = 10^2, lambda = 15), nrow = 10)
  splitting_results1 <- split_matrix(A, 'poisson', 0.3)
  splitting_results2 <- split_matrix(A, 'poisson', 0.5)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})
