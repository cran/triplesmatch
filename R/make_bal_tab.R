#' Make covariate balance
#'
#' @inheritParams boxplot_matches
#' @param X Covariate matrix
#' @param cov_names Row names to use instead of the column names of X when returning the table
#'
#' @return Table displaying means of the treated and control groups before and after matching,
#'     as well as standardized differences before and after matching
#' @export
#'
#' @examples
#' # Generate some data
#' set.seed(8)
#' n <- 200
#' nt <- floor(n * 0.5)
#' nc <- n - nt
#' x <- c(rnorm(nt, 0, 1), rnorm(nc, 0.6, 1))
#' z <- c(rep(1, nt), rep(0, nc))
#' names(z) <- 1:length(z)
#' names(x) <- 1:length(x)
#' # Create some strata
#' ps <- glm(z ~ x, family = binomial)$fitted.values
#' ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
#' # Create a distance matrix
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#' # Construct the triples match
#' triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
#' make_bal_tab(X = cbind(x, ps), z = z, m = triplesm$m, cov_names = c("x", "prop score"))

make_bal_tab <- function(X, z, m, cov_names) {
  if (is.null(names(z))) {
    stop("`z` must be a named vector")
  }
  if (is.null(row.names(X))) {
    stop("`X` must have row names")
  }
  mX <- makematrix(X, z, m)
  s1 <- apply(X[z == 1, ], 2, sd, na.rm=TRUE)
  s0 <- apply(X[z == 0, ], 2, sd, na.rm=TRUE)
  s <- sqrt(((s1^2) + (s0^2)) / 2)
  m1 <- apply(X[z == 1, ], 2, mean, na.rm=TRUE)
  m0 <- apply(X[z == 0, ], 2, mean, na.rm=TRUE)
  mean0N1 <- apply(mX[(mX$z == 0) & (mX$N == 1), 6:(dim(mX)[2])], 2, mean, na.rm=TRUE)
  mean0N2 <- apply(mX[(mX$z == 0) & (mX$N == 2), 6:(dim(mX)[2])], 2, mean, na.rm=TRUE)
  N1 <- sum(mX$N == 1)
  N2 <- sum(mX$N == 2)
  # Weighting the control means by the proportion of treated individuals in each type of triple
  m0w<-(N1 * mean0N1 + 2 * mean0N2 * N2) / (N1 + N2 * 2)
  o <- cbind(m1, m0, m0w, (m1-m0) / s, (m1-m0w) / s)
  colnames(o) <- c("Treated", "AllC", "MatchedC", "stdBefore", "stdAfter")
  rownames(o) <- c(cov_names)
  o
}
