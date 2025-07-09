#' In simulation settings where the true mean matrix is known and a set of
#' estimated communities are provided, this function returns the true target
#' of inference.
#'
#' @param M The true mean matrix.
#' @param u The linear combination vector (or matrix) which specifies which
#' connectivity parameters should be considered when constructing the selected
#' target of inference.
#' @param communities A vector or matrix which specifies the estimated
#' communities. If this is inputted as a vector, then if `M` is of size `n` x `n`,
#' then this should be a vector of length `n`, where the ith element is the
#' numbered community that the ith node belongs to. For example, `communities[2] = 3`
#' would indicate that the 2nd node belongs to the 3rd estimated community. If this
#' is inputted as a matrix, then it should be a matrix of 0s and 1s only, where
#' `communities[i, k] = 1` indicates that the ith node belongs to the kth
#' community. Each node is only allowed to belong to a single community, so
#' there should only be a single 1 in each row of this matrix.
#' @param K The number of estimated communities.
#' @param allow_self_loops A logical indicating whether the network allows
#' self loops (edges pointing from a node to itself.) By default this parameter
#' is set to `TRUE`. If this is set to `FALSE`, then the mean matrix should
#' have zeros along the diagonal.
#' @param is_directed A logical indicating whether the network is a directed
#' network, and by default is set to `TRUE`. If this is set to `FALSE`, then the
#' inputted mean matrix should be symmetric.
#' @param bernoulli_target A logical indicating whether we should calculate
#' the alternative target for Bernoulli networks specified in the paper that
#' is close to the actual target of interest. If this is set to `TRUE`, then
#' the argument `Atr` must also be provided.
#' @param gamma For Bernoulli networks, the parameter controlling the amount
#' of information allocated to the train network versus the test network. A
#' larger value of `gamma` indicates less information in the train network, and
#' more in the test network.
#' @param Atr A matrix of the same size as `M` arising from the
#' `split_network()` function. This argument will only be used if
#' the `bernoulli_target` argument is set to `TRUE`.
#' @returns A single numeric which is the true target of inference sought.
#' @examples
#' # As an example, let's build our matrix M to be a 10x10 matrix where each
#' # element is sampled uniformly from a Uniform(0, 1) distribution.
#' M <- matrix(stats::runif(n = 10^2), nrow = 10)
#'
#' # Set our "estimated_communities" to just have the first 5 nodes in the first
#' # community, and the last 5 nodes in the other community
#' communities_estimate <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
#'
#' # This choice of "u" vector would be conducting inference for the mean
#' # connectivity within the first estimated community.
#' u_vector <- c(1, 0,
#'               0, 0)
#'
#' # Conduct inference for the selected target (mean connectivity within the
#' # first estimated community)
#' target_of_inference <-
#'   check_target_of_inference(M = M, u = u_vector,
#'                             communities = communities_estimate, K = 2,
#'                             allow_self_loops = TRUE, is_directed = TRUE)
#' @export
check_target_of_inference <- function(M, u, communities, K,
                                      allow_self_loops = TRUE,
                                      is_directed = TRUE,
                                      bernoulli_target = FALSE,
                                      gamma = NULL,
                                      Atr = NULL) {

  # ---------------------- #
  # -- Check user input -- #
  # ---------------------- #

  if (!(is.matrix(M))) {
    stop("The input \"M\" needs to be a matrix.")
  }
  if (!(is.numeric(K))) {
    stop("The input \"K\" must be an integer.")
  }
  if (nrow(M) != ncol(M)) {
    stop("The input \"M\" must be a square matrix.")
  }
  if (is_directed == FALSE) {
    if (!isSymmetric(M)) {
      stop("When specifying an undirected network, the matrix \"M\" needs to be symmetric.")
    }
  }
  if (bernoulli_target) {
    if (is.null(Atr)) {
      stop("When setting \"bernoulli_target\" to TRUE, the argument \"Atr\" must be provided.")
    }
    if (is.null(gamma)) {
      stop("When setting \"bernoulli_target\" to TRUE, the argument \"gamma\" must be provided.")
    }
    if (!is.numeric(gamma)) {
      stop("The argument \"gamma\" must be numeric.")
    }
    if (!(is.matrix(Atr))) {
      stop("The input \"Atr\" needs to be a matrix.")
    }
    if (nrow(Atr) != ncol(Atr)) {
      stop("The input \"Atr\" must be a square matrix.")
    }
    if (nrow(Atr) != nrow(M)) {
      stop("The inputs \"M\" and \"Atr\" must be of the same dimensions.")
    }

  }

  # Extract some info from the matrices
  n <- NROW(M)
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
    # The number of nodes in each of the estimated communities
    n_hat <- apply(comm_as_matrix, 2, sum)
  }

  # =================================================== #
  # == Find out what the true target of inference is == #
  # =================================================== #

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
  if (!bernoulli_target) {
    # Where target is stored
    target_blocky <- matrix(rep(-999999, K^2), nrow = K)

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

          # Find the mean within this community pair
          mean_comm_pair <- 0
          for (ij_index in 1:NROW(Ikl)) {
            ij <- Ikl[ij_index, ]
            mean_comm_pair <- mean_comm_pair + M[ij[1], ij[2]]
          }
          target_blocky[k, l] <- mean_comm_pair / Ikl_card
        }
      }
    }

    # Based upon these results, return the final estimator as well as its variance
    target <- t(u) %*% as.vector(target_blocky)
    return(target)
  } else {

    # Convert the marginal mean matrix M into the conditional mean matrix T (Tmat)
    Cmask <- (gamma / (1-gamma))^(2*Atr - 1)
    Tmat <- M / (M + (1-M)*Cmask)

    # Calculate slightly shifted alternative target \xi from the paper
    target_blocky <- matrix(rep(-99999, K^2), nrow = K)

    for (k in 1:K) {
      for (l in 1:K) {
        # Only do calculation if corresponding entry of "u" is nonzero
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

          Ikl_card <- NROW(Ikl)

          # Calculate estimate
          # ------------------
          num_0s <- 0
          num_1s <- 0
          B0 <- 0
          B1 <- 0

          for (ij_index in 1:NROW(Ikl)) {
            ij <- Ikl[ij_index, ]

            if (Atr[ij[1], ij[2]] == 0) {
              B0 <- B0 + Tmat[ij[1], ij[2]]
              num_0s <- num_0s + 1
            } else if (Atr[ij[1], ij[2]] == 1) {
              B1 <- B1 + Tmat[ij[1], ij[2]]
              num_1s <- num_1s + 1
            } else {
              stop("Input \"Atr\" contains something other than a 0 or 1.")
            }
          }
          B0 <- B0 / num_0s
          B1 <- B1 / num_1s

          # Calculate Phi from B0 and B1
          Phi0 <- (num_0s / (num_0s + num_1s)) * (B0 / (B0 + ((1 - B0) * gamma / (1-gamma))))
          Phi1 <- (num_1s / (num_0s + num_1s)) * (B1 / (B1 + ((1 - B1) * (1-gamma) / gamma)))

          target_blocky[k, l] <- Phi0 + Phi1
        }
      }
    }

    # Based upon these results, return the final estimator as well as its variance
    target <- t(u) %*% as.vector(target_blocky)
    return(target)
  }
}
