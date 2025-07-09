#' Conduct inference using `Ate` for the selected target of inference: a
#' linear combination of connectivity parameters based upon estimated
#' communities. Note that communities should be estimated using the `Atr` matrix
#' produced from the `split_network()` function prior to using this function.
#'
#' @param Ate The test adjacency matrix produced from the `split_network()`
#' function which will be used to conduct inference.
#' @param u The linear combination vector (or matrix) which specifies which
#' connectivity parameters should be considered when constructing the selected
#' target of inference.
#' @param communities A vector or matrix which specifies the estimated
#' communities. If this is inputted as a vector, then if `Ate` is of size `n` x `n`,
#' then this should be a vector of length `n`, where the ith element is the
#' numbered community that the ith node belongs to. For example, `communities[2] = 3`
#' would indicate that the 2nd node belongs to the 3rd estimated community. If this
#' is inputted as a matrix, then it should be a matrix of 0s and 1s only, where
#' `communities[i, k] = 1` indicates that the ith node belongs to the kth
#' community. Each node is only allowed to belong to a single community, so
#' there should only be a single 1 in each row of this matrix.
#' @param distribution The distribution that the edges of the adjacency matrix
#' follow. Acceptable distributions are `"gaussian"`, `"poisson"`, or `"bernoulli"`.
#' @param K The number of estimated communities.
#' @param epsilon The parameter controlling the amount of information
#' allocated to the train network versus the test network. For Gaussian and
#' Poisson networks, a larger value of epsilon indicates more information in
#' the train network. For Bernoulli networks, this input is an alias to the
#' `gamma` parameter.
#' @param gamma For Bernoulli networks, the parameter controlling the amount
#' of information allocated to the train network versus the test network. A
#' larger value of `gamma` indicates less information in the train network, and
#' more in the test network.
#' @param Atr The train adjacency matrix produced from the `split_network()`
#' function. This is only necessary when the network has edges which follow
#' the Bernoulli distribution.
#' @param allow_self_loops A logical indicating whether the network allows
#' self loops (edges pointing from a node to itself.) By default this parameter
#' is set to `TRUE`. If this is set to `FALSE`, then the adjacency matrix should
#' have zeros along the diagonal.
#' @param is_directed A logical indicating whether the network is a directed
#' network, and by default is set to `TRUE`. If this is set to `FALSE`, then the
#' inputted adjacency matrix should be symmetric.
#' @param tau For networks with Gaussian edges only, this parameter indicates
#' the known common standard deviation (square root of the variance) of the
#' edges in the network.
#' @returns A list labeled with two elements labeled `"estimate"` and
#' `"estimate_variance"`, which contain the estimate for the selected target
#' of inference, as well as an estimate of the variance of the estimator.
#' @examples
#' # ==============================
#' # == Gaussian network example ==
#' # ==============================
#' # (Poisson networks proceed nearly identically except that the parameter
#' #  tau is not necessary to input into split_network() or infer_network())
#'
#' # First, split a simulated Gaussian adjacency matrix
#' A_gaussian <- matrix(stats::rnorm(n = 10^2, mean = 10, sd = 5), nrow = 10)
#' gaussian_split <- split_matrix(A = A_gaussian, distribution = "gaussian",
#'                                epsilon = 0.3, tau = 5)
#' A_gaussian_tr <- gaussian_split$Atr
#' A_gaussian_te <- gaussian_split$Ate
#'
#' # Estimate some communities using the train matrix using spectral
#' # clustering for K=3 communities
#' if (requireNamespace("nett", quietly = TRUE)) {
#'   communities_estimate <- nett::spec_clust(A_gaussian_tr, K = 3)
#' } else {
#'   communities_estimate <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
#' }
#'
#' # This particular "u" vector will specify that we want to conduct inference
#' # for the mean connectivity within the first estimated community.
#' # Note that "u" is of length 9, because we have 3 estimated communities.
#' # (This should always be of length K^3)
#' u_vector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'
#' # We can also specify "u" in matrix form.
#' u_matrix <- matrix(c(1, 0, 0,
#'                      0, 0, 0,
#'                      0, 0, 0), nrow = 3)
#'
#' # Conduct inference for the selected target (mean connectivity within the
#' # first estimated community)
#' gaussian_inference <-
#'     infer_network(Ate = A_gaussian_te, u = u_matrix,
#'                   communities = communities_estimate,
#'                   distribution = "gaussian",
#'                   epsilon = 0.3, K = 3, tau = 5)
#'
#' # Produce a 90% confidence interval for the target of inference
#' margin_of_error <- sqrt(gaussian_inference$estimate_variance) * qnorm(0.95)
#' ci_upper_bound <- gaussian_inference$estimate + margin_of_error
#' ci_lower_bound <- gaussian_inference$estimate - margin_of_error
#'
#' # ===============================
#' # == Bernoulli network example ==
#' # ===============================
#'
#' # First, split a simulated Bernoulli adjacency matrix with gamma=0.10
#' A_bernoulli <- matrix(stats::rbinom(n = 10^2, size = 1, p = 0.5), nrow = 10)
#' bernoulli_split <- split_matrix(A = A_bernoulli, distribution = "bernoulli",
#'                                 gamma = 0.10)
#' A_bernoulli_tr <- bernoulli_split$Atr
#' A_bernoulli_te <- bernoulli_split$Ate
#'
#' # Estimate some communities using the train matrix using spectral
#' # clustering for K=3 communities
#' if (requireNamespace("nett", quietly = TRUE)) {
#'   communities_estimate <- nett::spec_clust(A_bernoulli_tr, K = 3)
#' } else {
#'   communities_estimate <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
#' }
#'
#' # This particular "u" vector will specify that we want to conduct inference
#' # for the mean connectivity within the first estimated community.
#' # Note that "u" is of length 9, because we have 3 estimated communities.
#' # (This should always be of length K^3)
#' u_vector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'
#' # We can also specify "u" in matrix form.
#' u_matrix <- matrix(c(1, 0, 0,
#'                      0, 0, 0,
#'                      0, 0, 0), nrow = 3)
#'
#' # Conduct inference for the selected target (mean connectivity within the
#' # first estimated community)
#' bernoulli_inference <-
#'     infer_network(Ate = A_bernoulli_te, u = u_matrix,
#'                   communities = communities_estimate,
#'                   distribution = "bernoulli",
#'                   gamma = 0.10, K = 3,
#'                   Atr = A_bernoulli_tr)
#'
#' # Produce a 90% confidence interval for the target of inference
#' margin_of_error <- sqrt(bernoulli_inference$estimate_variance) * qnorm(0.95)
#' ci_upper_bound <- bernoulli_inference$estimate + margin_of_error
#' ci_lower_bound <- bernoulli_inference$estimate - margin_of_error
#' @export
infer_network <- function(Ate, u, communities, distribution, K, epsilon = 0.5,
                          gamma = NULL, Atr = NULL, allow_self_loops = TRUE, is_directed = TRUE,
                          tau = NA) {

  # ---------------------- #
  # -- Check user input -- #
  # ---------------------- #

  if (!(is.matrix(Ate))) {
    stop("The input \"Ate\" needs to be a matrix.")
  }
  if (!(is.numeric(K))) {
    stop("The input \"K\" must be an integer.")
  }
  if (nrow(Ate) != ncol(Ate)) {
    stop("The input \"Ate\" must be a square matrix.")
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
    if (is.na(tau)) {
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

    if (is.null(Atr)) {
      stop("When specifying a Bernoulli disribution, an input for \"Atr\" is required.")
    }
    if (!(is.matrix(Atr))) {
      stop("The input \"Atr\" into infer_matrix() needs to be a matrix.")
    }
    if (nrow(Atr) != ncol(Atr)) {
      stop("The input \"Atr\" must be a square matrix.")
    }
    if (any(!(Ate %in% c(0, 1)))) {
      stop("When specifying Bernoulli data, the matrix \"Ate\" must only contains 0s and 1s.")
    }
    if (any(!(Atr %in% c(0, 1)))) {
      stop("When specifying Bernoulli data, the matrix \"Atr\" must only contains 0s and 1s.")
    }
    if (is_directed == FALSE) {
      if (!isSymmetric(Atr)) {
        stop("When specifying an undirected network, the matrix \"Atr\" needs to be symmetric.")
      }
    }
  }
  if (is_directed == FALSE) {
    if (!isSymmetric(Ate)) {
      stop("When specifying an undirected network, the matrix \"Ate\" needs to be symmetric.")
    }
  }

  # Extract some info from the matrices
  n <- NROW(Ate)
  K <- as.integer(K) # Technically just needs to be input as a numeric

  # Check the communities input
  if (!((is.matrix(communities)) |
        (is.numeric(communities) && is.vector(communities)) |
        (is.factor(communities) && is.vector(communities)))) {
    stop("Input \"communities\" must be either a matrix or a vector.")
  }
  # If u is inputted as a matrix, first vectorize it
  if (is.matrix(u)) {
    u_as_matrix <- u
    u <- as.vector(u)
  } else {
    u_as_matrix <- matrix(u, ncol = K)

    # Check the linear combination vector u
    if (!(is.vector(u))) {
      stop("The linear combination vector \"u\" must be a vector.")
    }
    if (!(is.numeric(u))) {
      stop("The linear combination vector \"u\" must be numeric.")
    }
    if (length(u) != K^2) {
      stop("The length of the linear combination vector \"u\" must be K^2.")
    }
  }
  # Make sure that "u" is of norm 1
  u <- u / sqrt(sum(u^2))
  # If communities is a vector, then convert it to a matrix as well.
  if (is.vector(communities)) {
    if (!is.factor(communities) & !is.numeric(communities)) {
      stop("When inputting \"communities\" as a vector, it must be a vector of factors or a vector of integers.")
    }
    if (length(communities) != n) {
      stop("The length of the vector \"communities\" is not equal to the dimension of the network \"Ate\".")
    }
    if (is.numeric(communities)) {
      communities <- as.integer(communities)
    }

    comm_as_vector <- communities
    comm_as_matrix <- matrix(rep(NA, n*K), nrow = n)
    for (i in 1:K) {
      comm_as_matrix[, i] <- 1 * (comm_as_vector == i)
    }
    n_hat <- apply(comm_as_matrix, 2, sum) # The number of nodes in each of the
                                           # estimated communities
  }
  if (is.matrix(communities)) {
    # Check right dimensions
    if (NROW(communities) != n) {
      stop("The number of rows in \"communities\" matrix must be equal to the nodes in the network.")
    }
    if (NCOL(communities) != K) {
      stop("The number of columns in \"communities\" matrix must be equal to the number of estimated communities \"K\".")
    }
    if (any(!(communities %in% c(0, 1)))) {
      stop("When specifying \"communities\" as a matrix, it must only contain 0s and 1s.")
    }
    if (any(apply(communities, 1, sum) != 1)) {
      stop("When specifying \"communities\" as a matrix, each row must only contain a single 1. No partial community membership is allowed.")
    }
    comm_as_matrix <- communities
    comm_as_vector <- rep(NA, n)
    for (i in 1:n) {
      comm_as_vector[i] <- which(communities[i, ] == 1)
    }
    n_hat <- apply(comm_as_matrix, 2, sum) # The number of nodes in each of the
                                           # estimated communities
  }

  # ======================================================= #
  # == Do inference depending on what distribution it is == #
  # ======================================================= #

  # -------------------------------------------------------------------------- #
  # -- Some notes on implementation: there are some ways to make the        -- #
  # -- code really elegant with nice linear algebra expressions that run    -- #
  # -- quickly. However, they don't really work for undirected networks or  -- #
  # -- networks without self loops. So the choice is between writing a huge -- #
  # -- amount of code that does a case-by-case basis but runs faster, or    -- #
  # -- writing code that works for all situations but is slower. I chose to -- #
  # -- do the latter.                                                       -- #
  # -------------------------------------------------------------------------- #

  # This code is fundamentally similar whether it is Gaussian, Poisson, or
  # Bernoulli edges.

  # Where estimator and variance is stored
  estimator_blocky <- matrix(rep(-999999, K^2), nrow = K)
  estimator_variance_blocky <- matrix(rep(-999999, K^2), nrow = K)

  # Repeat this process for all the community pairs, but only the ones
  # that correspond to a nonzero entry for u.
  for (k in 1:K) {
    for (l in 1:K) {
      # Only do the calculation if the corresponding entry of u
      # is nonzero.
      if (u_as_matrix[k, l] != 0) {
        k_nodes <- which(comm_as_vector == k)
        l_nodes <- which(comm_as_vector == l)

        # Refine our index set
        if (is_directed) {
          Ikl <- as.matrix(expand.grid(k_nodes, l_nodes))
        } else {
          Ikl <- rbind(as.matrix(expand.grid(k_nodes, l_nodes)),
                       as.matrix(expand.grid(l_nodes, k_nodes)))
          Ikl <- Ikl[Ikl[, 1] <= Ikl[, 2], ]
        }
        if (!allow_self_loops) {
          Ikl <- Ikl[Ikl[, 1] != Ikl[, 2], ]
        }

        # Count up cardinality
        Ikl_card <- NROW(Ikl)

        # -------------- #
        # -- Gaussian -- #
        # -------------- #

        if (distribution == "gaussian") {
          estimator_blocky_k_l <- 0
          for (ij_index in 1:NROW(Ikl)) {
            ij <- Ikl[ij_index, ]
            estimator_blocky_k_l <- estimator_blocky_k_l + Ate[ij[1], ij[2]]
          }
          estimator_blocky[k, l] <- (estimator_blocky_k_l / Ikl_card) / (1-epsilon)
          estimator_variance_blocky[k, l] <- (tau^2 / Ikl_card) / (1-epsilon)
        }

        # ------------- #
        # -- Poisson -- #
        # ------------- #

        if (distribution == "poisson") {
          estimator_blocky_k_l <- 0
          for (ij_index in 1:NROW(Ikl)) {
            ij <- Ikl[ij_index, ]
            estimator_blocky_k_l <- estimator_blocky_k_l + Ate[ij[1], ij[2]]
          }
          estimator_blocky[k, l] <- (estimator_blocky_k_l / Ikl_card) / (1-epsilon)
          estimator_variance_blocky[k, l] <- (estimator_blocky[k, l] / Ikl_card) / (1-epsilon)
        }

        # --------------- #
        # -- Bernoulli -- #
        # --------------- #

        if (distribution == "bernoulli") {

          # Some constants
          gamma <- epsilon # Just because this is the notation used in the paper
          c0 <- log(gamma / (1-gamma))
          c1 <- log((1-gamma) / gamma)

          B0 <- 0
          B1 <- 0
          num_0s <- 0
          num_1s <- 0

          for (ij_index in 1:NROW(Ikl)) {
            ij <- Ikl[ij_index, ]

            if (Atr[ij[1], ij[2]] == 0) {
              B0 <- B0 + Ate[ij[1], ij[2]]
              num_0s <- num_0s + 1
            } else {
              B1 <- B1 + Ate[ij[1], ij[2]]
              num_1s <- num_1s + 1
            }
          }
          B0 <- B0 / num_0s
          B1 <- B1 / num_1s

          Phi_kl_0 <- 0
          Phi_kl_1 <- 0
          Delta_kl_0 <- 0
          Delta_kl_1 <- 0

          if (num_0s > 0) {
            # Mean estimate
            Phi_kl_0 <- (num_0s / Ikl_card) * (B0 / (B0 + ((1 - B0) * gamma / (1-gamma))))

            # Variance estimate
            B0_kl_adjust <- B0
            if (B0 == 0) {
              B0_kl_adjust <- 1 / (2*num_0s)
            }
            if (B0 == 1) {
              B0_kl_adjust <- (num_0s - 0.5) / num_0s
            }
            Delta_kl_0 <- (B0_kl_adjust * (1 - B0_kl_adjust) * exp(2*c0)) /
              (num_0s * ((1 - B0)*exp(c0) + B0)^4)
          }
          if (num_1s > 0) {
            # Mean estimate
            Phi_kl_1 <- (num_1s / Ikl_card) * (B1 / (B1 + ((1 - B1) * (1-gamma) / gamma)))

            # Variance estimate
            B1_kl_adjust <- B1
            if (B1 == 0) {
              B1_kl_adjust <- 1 / (2*num_1s)
            }
            if (B1 == 1) {
              B1_kl_adjust <- (num_1s - 0.5) / num_1s
            }
            Delta_kl_1 <- (B1_kl_adjust * (1 - B1_kl_adjust) * exp(2*c1)) /
              (num_1s * ((1 - B1)*exp(c1) + B1)^4)
          }

          # Bound these variance estimates by 0.25 if necessary
          Delta_kl_0 <- pmin(Delta_kl_0, 0.25)
          Delta_kl_1 <- pmin(Delta_kl_1, 0.25)

          # Store as final mean and variance estimate
          estimator_blocky[k, l] <- Phi_kl_0 + Phi_kl_1
          estimator_variance_blocky[k, l] <- (num_0s / Ikl_card)^2 * Delta_kl_0 +
            (num_1s / Ikl_card)^2 * Delta_kl_1
        }
      }
    }
  }

  # Based upon these results, return the final estimator as well as its variance
  estimator <- t(u) %*% as.vector(estimator_blocky)
  estimator_variance <- t(u) %*% diag(as.vector(estimator_variance_blocky)) %*% u
  return(list("estimate" = estimator, "estimate_variance" = estimator_variance))
}
