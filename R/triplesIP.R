#' Integer program for finding optimal triples match
#'
#' This finds the optimal triples match using a quadratic program. The 'gurobi'
#' package should be installed if using this function. This function should not be
#' used for large problems. Note that this solver may find a good solution even if not
#' optimal; setting `time_limit` is recommended. For most problems, [triples()] should
#' be used instead to find a good approximate solution very quickly.
#'
#' @param z Treatment indicator vector. 0 for control, 1 for treated
#' @param cost Matrix of costs. Rows correspond to treated units; columns to controls
#' @param mt The number of treated units to be used
#' @param mc The number of control units to be used
#' @param time_limit The amount of time in seconds before the solver should abort
#' @param threads The number of threads that should be allocated
#' @param verbose Whether the output of the 'gurobi' solver should be printed. 0 if not, 1 if so
#'
#' @return A named list with two elements: `match` and `opt_info`. `match` contains the triples match.
#' Similarly to the [triples()] function, this is in the form of a
#'  data.frame with number of rows equal to the number of triples and 8 columns specifying the
#'  match number, the names of the three units within the match, the costs of the two
#'  treated-control pairs within the match, the number of treated units, and the stratum.
#'  `opt_info` contains technical output from the optimization solver.
#' @export
#' @importFrom rlang is_installed
#'
#' @seealso triples for an approximate solution
#'
#' @examplesIf rlang::is_installed(c("gurobi"))
#' # Generate some data
#' set.seed(1)
#' n <- 40
#' x <- rnorm(n, 0, 1)
#' nt <- floor(n * 0.4)
#' nc <- n - nt
#' z <- c(rep(1, nt), rep(0, nc))
#' # Create some strata
#' ps <- glm(z ~ x, family = binomial)$fitted.values
#' ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
#' # Create a distance matrix
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#' # Construct the triples match using integer program for stratum 1
#' mIP <- triplesIP(z = z[ps_st == 1], cost = dist[[1]],
#'                  mt = 5, mc = 7, time_limit = 30, threads = 1, verbose = 0)
#'

triplesIP <- function(z, cost, mt, mc, time_limit = Inf, threads = 1, verbose = 0) {

  if(!requireNamespace("gurobi", quietly = TRUE)) {
    stop("'gurobi' must be installed to use this function. See triples() function for an alternative. \n")
  }

  nt <- sum(z)
  n <- length(z)
  nc <- n - nt
  stopifnot(nt == nrow(cost))
  stopifnot(nc == ncol(cost))

  ntrip <- (mt+mc) / 3

  if(nt == 0) {
    warning("There are no treated individuals to be matched in this stratum.")
    return(NULL)
  }
  if(nc == 0) {
    warning("There are no controls to be matched in this stratum.")
    return(NULL)
  }
  if(ntrip == 0) {
    warning("There are not enough individuals to match in triples.")
    return(NULL)
  }


  # Each person has an indicator for each triple,
  #     and then there is a slack variable for each triple
  nvar <- n * ntrip + ntrip

  Q <- matrix(0, nrow = nvar, ncol = nvar)
  # For each triple, the quadratic matrix in the objective
  #    is the cost matrix in the upper right corner of a nxn matrix of 0s
  for (i in 1:ntrip) {
    Q[(i-1)*n + 1:nt, (i-1)*n + ((nt+1):n)] <- cost
  }

  # There are constraints for the total number of units in each triple,
  #     the number of treated in each triple,
  #     the total number of treated,
  #     the total number of control,
  #     and the number of times each unit is used
  ncon <- 2 * ntrip + 2 + n
  b <- rep(0, ncon)
  A <- matrix(0, nrow = ncon, ncol = nvar)
  for (i in 1:ntrip) {
    # Total number of units in each triple should be 3
    A[i, (i-1)*n + 1:n] <- 1
    b[i] <- 3
    # Number of treated plus a slack variable in each triple should be 2
    A[ntrip + i, c((i-1)*n + 1:nt, (n * ntrip + i))] <- 1
    b[ntrip + i] <- 2
    # Sum the number of treated across triples
    A[2 * ntrip + 1, (i-1)*n + 1:nt] <- 1
    # Sum the number of controls across triples
    A[2 * ntrip + 2, (i-1)*n + (nt+1):n] <- 1
  }
  b[2 * ntrip + 1] <- mt
  b[2 * ntrip + 2] <- mc
  for (i in 1:n) {
    A[2 * ntrip + 2 + i, (0:(ntrip - 1)) * n + i] <- 1
  }
  b[2 * ntrip + 2 + 1:n] <- 1

  lb <- rep(0, nvar)
  ub <- rep(1, nvar)

  int_con <- rep(TRUE, nvar)

  # All linear constraints are equalities except the final n which are \leq
  # Gurobi uses < to mean \leq
  sense <- c(rep("=", ncon - n), rep("<", n))

  model <- list(Q = Q,    # Quadratic objective matrix
                A = A,    # Linear constraint matrix
                sense = sense,    # Type of linear inequality
                rhs = b,  # Right hand side of linear constraints
                lb = lb,  # Lower bound of variables
                ub = ub,  # Upper bound of variables
                vtype = "I"   # All variables are integer
                )
  o <- gurobi::gurobi(model, params = list(TimeLimit = time_limit, Threads = threads, OutputFlag = verbose, MIPFocus = 0))

  # Format the results
  res <- data.frame(matrix(NA, nrow = ntrip, ncol = 7,
                           dimnames = list(NULL, c("triple", "itreated", "jcontrol", "kthird",
                                      "costStep1", "costStep2", "nOfTreated"))))

  for (i in 1:ntrip) {
    res[i, "triple"] <- i
    t_used <- which(o$x[(i-1)*n + 1:nt] == 1)
    c_used <- which(o$x[(i-1)*n + nt + 1:nc] == 1)
    res[i, "itreated"] <- rownames(cost)[t_used[1]]
    res[i, "jcontrol"] <- colnames(cost)[c_used[1]]
    res[i, "costStep1"] <- cost[res[i, "itreated"], res[i, "jcontrol"]]
    if (length(t_used) > 1) {
      res[i, "kthird"] <- rownames(cost)[t_used[2]]
      res[i, "nOfTreated"] <- 2
      res[i, "costStep2"] <- cost[res[i, "kthird"], res[i, "jcontrol"]]
    } else {
      res[i, "kthird"] <- colnames(cost)[c_used[2]]
      res[i, "nOfTreated"] <- 1
      res[i, "costStep2"] <- cost[res[i, "itreated"], res[i, "kthird"]]
    }
  }

  return(list(match = res, opt_info = o))
}
