#' Split an adjacency matrix into a train and test adjacency matrix using
#' either data thinning (Poisson or Gaussian edges) or data fission (Bernoulli
#' edges.)
#'
#' @param A The square adjacency matrix to be split.
#' @param distribution The distribution that the edges of the adjacency matrix
#' follow. Acceptable distributions are `"gaussian"`, `"poisson"`, or `"bernoulli"`.
#' @param epsilon The parameter controlling the amount of information
#' allocated to the train network versus the test network. For Gaussian and
#' Poisson networks, a larger value of `epsilon` indicates more information in
#' the train network. For Bernoulli networks, this input is an alias to the
#' `gamma` parameter.
#' @param gamma For Bernoulli networks, the parameter controlling the amount
#' of information allocated to the train network versus the test network. A
#' larger value of `gamma` indicates less information in the train network, and
#' more in the test network.
#' @param allow_self_loops A logical indicating whether the network allows
#' self loops (edges pointing from a node to itself.) By default this parameter
#' is set to `TRUE`. If this is set to `FALSE`, then the values in the
#' adjacency matrix along the diagonal will be ignored.
#' @param is_directed A logical indicating whether the network is a directed
#' network, and by default is set to `TRUE`. If this is set to `FALSE`, then
#' only the values along the upper triangular portion of the matrix will be used.
#' @param tau For networks with Gaussian edges only, this parameter indicates
#' the known common standard deviation (square root of the variance) of the
#' edges in the network.
#' @returns A list labeled with two elements labeled `"Atr"` and `"Ate"`, which
#' are the train and test networks, respectively.
#' @examples
#' # Split a simulated Gaussian adjacency matrix
#' A_gaussian <- matrix(rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
#' gaussian_split <- split_matrix(A_gaussian, "gaussian", 0.3, tau = 5)
#' A_gaussian_tr <- gaussian_split$Atr
#' A_gaussian_te <- gaussian_split$Ate
#'
#' # Split a simulated Bernoulli adjacency matrix with gamma = 0.25
#' A_bernoulli <- matrix(rbinom(n = 10^2, size = 1, p = 0.5), nrow = 10)
#' bernoulli_split <- split_matrix(A_bernoulli, "bernoulli", gamma = 0.25)
#' A_bernoulli_tr <- bernoulli_split$Atr
#' A_bernoulli_te <- bernoulli_split$Ate
#' @export
split_matrix <- function(A, distribution, epsilon = 0.5, gamma = NULL,
                         allow_self_loops = TRUE, is_directed = TRUE,
                         tau = NULL) {

  # ---------------------- #
  # -- Check user input -- #
  # ---------------------- #

  if (!(is.matrix(A))) {
    stop("The input \"A\" into split_matrix() needs to be a matrix.")
  }
  if (nrow(A) != ncol(A)) {
    stop("The input A must be a square matrix.")
  }
  if (!(distribution %in% c("gaussian", "poisson", "bernoulli"))) {
    stop("Input \"distribution\" needs to be one of \"gaussian\", \"poisson\", or \"bernoulli\".")
  }
  if (distribution %in% c("gaussian", "poisson")) {
    if (!is.null(gamma)) {
      warning("A value of \"gamma\" was input even though the distribution is Poisson or Gaussian. The value of \"gamma\" will be ignored.")
    }
  }
  if (distribution == "gaussian") {
    if (is.null(tau)) {
      stop("When specifying a Gaussian distribution, the known standard deviation \"tau\" must be inputted.")
    }
    if (!is.numeric(tau)) {
      stop("When specifying a Gaussian distribution, the inputted standard deviation \"tau\" needs to be numeric and positive.")
    }
    if (tau <= 0) {
      stop("When specifying a Gaussian distribution, the inputted \"tau\" needs to be positive.")
    }
  }
  if (distribution == "bernoulli") {
    # Assign the input gamma to epsilon if it exists
    if (!is.null(gamma)) {
      if (!is.numeric(gamma)) {
        stop("The input \"gamma\" must be numeric.")
      }
      epsilon <- gamma
    }
  }

  # ---------------------------------------------------------
  # -- Do thinning / fission depending on the distribution --
  # ---------------------------------------------------------

  n <- nrow(A)

  if (distribution == "poisson") {
    A_tr <- matrix(stats::rbinom(n = n^2, size = as.vector(A), prob = epsilon), nrow = n)
    A_te <- A - A_tr
  }
  if (distribution == "gaussian") {
    A_tr <- matrix(stats::rnorm(n = n^2, mean = epsilon * as.vector(A),
                         sd = sqrt(epsilon*(1-epsilon)) * tau), nrow = n)
    A_te <- A - A_tr
  }
  if (distribution == "bernoulli") {
    Zfission <- matrix(stats::rbinom(n = n^2, size = 1, prob = epsilon), nrow = n)
    A_tr <- A * (1 - Zfission) + (1 - A) * Zfission
    A_te <- A
  }

  # After doing the above, if we have an undirected matrix, then clear out the
  # lower triangular portion and replace it with the upper triangular portion.
  # Similarly, if the matrix does not allow self-loops then clear it out and
  # replace it with zeros.
  if (!is_directed) {
    A_tr[lower.tri(A_tr)] <- 0
    A_te[lower.tri(A_te)] <- 0
  }
  if (!allow_self_loops) {
    diag(A_tr) <- 0
    diag(A_te) <- 0
  }

  return(list("Atr" = A_tr, "Ate" = A_te))
}
