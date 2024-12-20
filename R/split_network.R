split_matrix <- function(A, distribution, epsilon, allow_self_loops = TRUE,
                         is_directed = TRUE, sigma = NA) {

  # ----------------------
  # -- Check user input --
  # ----------------------

  if (!('matrix' %in% class(A))) {
    stop('The input A into split_matrix() needs to be a matrix.')
  }
  if (nrow(A) != ncol(A)) {
    stop('The input A must be a square matrix.')
  }
  if (!(distribution %in% c('gaussian', 'poisson', 'bernoulli'))) {
    stop('Input \'distribution\' needs to be one of \'gaussian\', \'poisson\', or \'bernoulli\'.')
  }
  if (distribution == 'gaussian') {
    if (is.na(sigma)) {
      stop('When specifying a Gaussian distribution, a \'sigma\' must be inputted.')
    }
    if (!is.numeric(sigma)) {
      stop('When specifying a Gaussian distribution, the inputted sigma needs to be numeric and positive.')
    }
    if (sigma <= 0) {
      stop('When specifying a Gaussian distribution, the inputted sigma needs to be positive.')
    }
  }
  if (is_directed == FALSE) {
    if (!isSymmetric(A)) {
      stop('When specifying an undirected network, the matrix A needs to be symmetric.')
    }
  }

  # ---------------------------------------------------------
  # -- Do thinning / fission depending on the distribution --
  # ---------------------------------------------------------

  n <- nrow(A)

  if (distribution == 'poisson') {
    A_tr <- matrix(rbinom(n = n^2, size = as.vector(A), prob = epsilon), nrow = n)
    A_te <- A - A_tr
  }
  if (distribution == 'gaussian') {
    A_tr <- matrix(rnorm(n = n^2, mean = epsilon * as.vector(A),
                         sd = sqrt(epsilon*(1-epsilon)) * sigma), nrow = n)
    A_te <- A - A_tr
  }
  if (distribution == 'bernoulli') {
    Zfission <- matrix(rbinom(n = n^2, size = 1, prob = epsilon), nrow = n)
    A_tr <- A * (1 - Zfission) + (1 - A) * Zfission
    A_te <- A
  }

  # After doing the above, if we have an undirected matrix, then clear out the
  # lower triangular portion and replace it with the upper triangular portion.
  # Similarly, if the matrix does not allow self-loops then clear it out and
  # replace it with zeros.
  if (!is_directed) {
    A_tr[lower.tri(A_tr)] <- t(A_tr)[lower.tri(A_tr)]
    A_te[lower.tri(A_te)] <- t(A_te)[lower.tri(A_te)]
  }
  if (!allow_self_loops) {
    diag(A_tr) <- 0
    diag(A_te) <- 0
  }

  return(list('Atr' = A_tr, 'Ate' = A_te))
}
