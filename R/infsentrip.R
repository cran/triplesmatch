#' Sensitivity analysis for triples matches informed by tests for unmeasured bias
#'
#' This function is very similar to [informedSen::informedsen()], with a few minor differences.
#' This version allows for matches to contain either two treated units and one control unit
#' or two controls and one treated unit
#' ([informedSen::informedsen()] requires only one treated unit and a fixed number of controls in each match).
#' This version also allows the primary hypothesis test to be one-sided.
#' To use this function, the optimization software 'gurobi' and its R package must be installed.
#'
#' @inheritParams sentrip
#' @param sc A matrix with N rows and at least two columns. The first column is the primary outcome
#' and the remaining columns are unaffected outcomes used to test for bias
#' @param z A vector of length N with 0s for control and 1s for treated units
#' @param mset A vector of length N indicating the matched triple
#' @param alpha A vector with length equal to the number of columns of `sc`. Each coordinate
#' contains the level of the test applied to the corresponding column of `sc`. If `alpha`
#' is a scalar, it is repeated for each column
#' @param alternative One of `greater`, `less` or `both`. `greater` implies the
#' alternative hypothesis is that the treatment has a positive effect on the scores of the primary outcome,
#' `less` implies the alternative hypothesis is that the treatment has a negative
#' effect on the scores of the primary outcome, and `both` conducts a two-sided hypothesis test.
#' The negative outcomes will always be two-sided tests (since one does not expect an effect in either direction)
#'
#' @return \describe{
#' \item{result }{Text indicating whether or not the test for bias rejects all biases of magnitude Gamma or less.  If yes, then the conclusion is that you must increase Gamma to continue.  If no, then the test on the primary outcome is conducted inside the confidence set defined by a test for bias.}
#' \item{optimization.problem }{Reiterates the result above, where the word yes means the optimization problem is infeasible, and the word no means it is feasible.  See the conclusion for a scientific interpretation of this aspect of the output.}
#' \item{conclusion }{Text indicating the result of the test for effect on the primary outcome.}
#' \item{deviates }{A vector of standardized deviates that might be compared with the standard Normal distribution.  There is one deviate for each column of `sc`.  If `sc` has column names, then the column names label the
#'   deviates.  The deviates are computed at the treatment assignment probabilities, theta, that solve the constrained optimization problem.}
#' \item{alphas }{A vector of significance levels used for the deviates, together with their total.  The total is relevant if the Bonferroni inequality is used to ensure joint level of all the tests.  The absolute deviates might be compared with qnorm(1-alphas/2) for a two-sided test.}
#' }
#' @seealso formattrip for easily creating inputs to this function.
#' @export
#' @importFrom rlang is_installed
#'
#' @examplesIf rlang::is_installed(c("gurobi", "Matrix"))
#' # Generate some data
#' set.seed(316)
#' n <- 30
#' x <- rnorm(n, 0, 1)
#' nt <- floor(n * 0.2)
#' nc <- n - nt
#' z <- c(rep(1, nt), rep(0, nc))
#' # Create a distance matrix (all units in one stratum here)
#' dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
#' # Conduct the triples match
#' triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")
#' # Create primary and negative outcomes with some random unit names
#' y <- cbind(rnorm(40), runif(40))
#' rownames(y) <- sample(1:40)
#' # Reformat the triples match
#' ylong <- formattrip(m = triplesm, y = y, type = "long")
#' # Score the outcomes, in this case aberrant ranks for primary outcome and
#' #      ranks for unaffected outcome
#' y[, 1] <- aberrantscoreslong(y[, 1], cutoff = 0.5, cutoff_dir = "greater")
#' y[, 2] <- rank(y[, 2])
#' # Run the informed sensitivity analysis at gamma of 1.5
#' inf1 <- infsentrip(gamma = 1.5, sc = ylong$y, z = ylong$z, ylong$mset,
#'                   alpha = 0.05, alternative = "both")

infsentrip <- function(gamma, sc, z, mset, alpha, alternative = "both") {
  if(!requireNamespace("gurobi", quietly = TRUE)) {
    stop("gurobi must be installed to use this function. \n")
  }
  if(!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix must be installed to use this function. \n")
  }

  if (length(alpha) == 1) {
    alpha <- rep(alpha, ncol(sc))
  }
  stopifnot(length(alpha) == ncol(sc))
  stopifnot(length(z) == nrow(sc))
  stopifnot(length(mset) == nrow(sc))
  stopifnot(gamma >= 1)
  stopifnot(alternative %in% c("less", "greater", "both"))

  rownames(sc) <- NULL

  # Convert manyT:1C to 1T:manyC
  # Multiply scores by -1 and change the roles
  for (i in unique(mset)) {
    if (sum(z[mset == i]) > 1) {
      sc[mset == i, ] <- -sc[mset == i, ]
      match_total <- colSums(sc[mset == i, ])
      t_ind <- which(z == 1 & mset == i)
      c_ind <- which(z == 0 & mset == i)
      z[c_ind] <- 1
      z[t_ind] <- 0
    }
  }

  # For the primary outcome, to have a one sided test at alpha,
  #    we input alpha*2 for the two sided test.
  #    If the test rejects, we will check whether the deviate is in the correct one-sided direction
  #    before confirming the results.
  if (alternative != "both") {
    alpha[1] <- alpha[1] * 2
  }

  invisible(utils::capture.output({
    infsen <- informedSen::informedsen(gamma = gamma, sc = sc / 100,
                                       z = z, mset = mset,
                                       alpha = alpha)
  }))

  # Change the alpha level in the conclusion back to the one sided level
  if (alternative != "both") {
    infsen$conclusion <- gsub(alpha[1], alpha[1]/2, infsen$conclusion)
    infsen$alphas[1] <- infsen$alpha[1]/2
    infsen$alphas['total'] <- sum(infsen$alphas[-length(infsen$alphas)])
  }

  # If the test rejects but the deviate is in the wrong direction, change the conclusion
  if (alternative == "greater" & !is.na(infsen$deviates[1]) & infsen$deviates[1] < 0) {
    infsen$conclusion <- gsub("rejects", "fails to reject", infsen$conclusion)
  } else if (alternative == "less" & !is.na(infsen$deviates[1]) & infsen$deviates[1] > 0) {
    infsen$conclusion <- gsub("rejects", "fails to reject", infsen$conclusion)
  }

  return(infsen)
}
