# -------------------- #
# -- Error handling -- #
# -------------------- #

test_that('Only accepts a matrix as input', {
  A <- "I like ice cream"
  expect_error(split_network(A, 'gaussian', 0.3, tau = 5))
})

test_that('Only accepts square matrices', {
  A <- matrix(rnorm(n = 4*3, mean = 0, sd = 1), nrow = 4, ncol = 3)
  expect_error(split_network(A, 'gaussian', 0.3, tau = 5))
})

test_that('Only accepts the right type of distributions', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  expect_error(split_network(A, 'dirichlet', 0.3, tau = 5))
})

# -------------------- #
# -- Correct output -- #
# -------------------- #

# Gaussian tests
test_that('Gaussian splitting works (directed network with self loops)', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  splitting_results1 <- split_network(A, 'gaussian', 0.3, tau = 5)
  splitting_results2 <- split_network(A, 'gaussian', 0.5, tau = 5)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})

test_that('Gaussian splitting works (directed network without self loops)', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  diag(A) <- 0
  splitting_results1 <- split_network(A, 'gaussian', 0.3,
                                     allow_self_loops = FALSE, tau = 5)
  splitting_results2 <- split_network(A, 'gaussian', 0.5,
                                     allow_self_loops = FALSE, tau = 5)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})

test_that('Gaussian splitting works (undirected network with self loops)', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  A <- A + t(A)
  splitting_results1 <- split_network(A, 'gaussian', 0.3,
                                     is_directed = FALSE, tau = 5)
  splitting_results2 <- split_network(A, 'gaussian', 0.5,
                                     is_directed = FALSE, tau = 5)

  recovered_matrix_1 <- splitting_results1$Atr + splitting_results1$Ate
  recovered_matrix_2 <- splitting_results2$Atr + splitting_results2$Ate

  # Check that we recover the original matrix
  expect_equal(A[upper.tri(A)], recovered_matrix_1[upper.tri(recovered_matrix_1)])
  expect_equal(A[upper.tri(A)], recovered_matrix_2[upper.tri(recovered_matrix_2)])
})

test_that('Gaussian splitting works (undirected network without self loops)', {
  A <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
  A <- A + t(A)
  diag(A) <- 0
  splitting_results1 <- split_network(A, 'gaussian', 0.3,
                                     allow_self_loops = FALSE,
                                     is_directed = FALSE, tau = 5)
  splitting_results2 <- split_network(A, 'gaussian', 0.5,
                                     allow_self_loops = FALSE,
                                     is_directed = FALSE, tau = 5)

  recovered_matrix_1 <- splitting_results1$Atr + splitting_results1$Ate
  recovered_matrix_2 <- splitting_results2$Atr + splitting_results2$Ate

  # Check that we recover the original matrix
  expect_equal(A[upper.tri(A)], recovered_matrix_1[upper.tri(recovered_matrix_1)])
  expect_equal(A[upper.tri(A)], recovered_matrix_2[upper.tri(recovered_matrix_2)])
})

# Poisson tests
test_that('Poisson splitting works (directed network with self loops)', {
  A <- matrix(rpois(n = 10^2, lambda = 15), nrow = 10)
  splitting_results1 <- split_network(A, 'poisson', 0.3)
  splitting_results2 <- split_network(A, 'poisson', 0.5)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})

test_that('Poisson splitting works (directed network without self loops)', {
  A <- matrix(rpois(n = 10^2, lambda = 15), nrow = 10)
  diag(A) <- 0
  splitting_results1 <- split_network(A, 'poisson', 0.3,
                                     allow_self_loops = FALSE)
  splitting_results2 <- split_network(A, 'poisson', 0.5,
                                     allow_self_loops = FALSE)

  expect_equal(A, splitting_results1$Atr + splitting_results1$Ate)
  expect_equal(A, splitting_results2$Atr + splitting_results2$Ate)
})

test_that('Poisson splitting works (undirected network with self loops)', {
  A <- matrix(rpois(n = 10^2, lambda = 15), nrow = 10)
  A <- A + t(A)
  splitting_results1 <- split_network(A, 'poisson', 0.3,
                                     is_directed = FALSE)
  splitting_results2 <- split_network(A, 'poisson', 0.5,
                                     is_directed = FALSE)

  recovered_matrix_1 <- splitting_results1$Atr + splitting_results1$Ate
  recovered_matrix_2 <- splitting_results2$Atr + splitting_results2$Ate

  # Check that we recover the original matrix
  expect_equal(A[upper.tri(A)], recovered_matrix_1[upper.tri(recovered_matrix_1)])
  expect_equal(A[upper.tri(A)], recovered_matrix_2[upper.tri(recovered_matrix_2)])
})

test_that('Poisson splitting works (undirected network without self loops)', {
  A <- matrix(rpois(n = 10^2, lambda = 15), nrow = 10)
  A <- A + t(A)
  diag(A) <- 0
  splitting_results1 <- split_network(A, 'poisson', 0.3,
                                     allow_self_loops = FALSE,
                                     is_directed = FALSE)
  splitting_results2 <- split_network(A, 'poisson', 0.5,
                                     allow_self_loops = FALSE,
                                     is_directed = FALSE)

  recovered_matrix_1 <- splitting_results1$Atr + splitting_results1$Ate
  recovered_matrix_2 <- splitting_results2$Atr + splitting_results2$Ate

  # Check that we recover the original matrix
  expect_equal(A[upper.tri(A)], recovered_matrix_1[upper.tri(recovered_matrix_1)])
  expect_equal(A[upper.tri(A)], recovered_matrix_2[upper.tri(recovered_matrix_2)])
})
