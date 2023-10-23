#' Make Mahalanobis distance matrix
#'
#' @param X Covariate matrix to be used in calculating distances
#' @param z Vector of treatment assignments
#' @param rank_based Whether to use rank based Mahalanobis distance or not
#' @param st Vector of stratum assignments, should be denoted by consecutive integers starting from 1
#'
#' @return List of squared Mahalanobis distance matrices between each pair of treated-control units in a stratum.
#'      There is one entry in the list for each stratum.
#'      This entry is a distance matrix with a row for each treated unit and a column
#'      for each control unit in the stratum.
#' @export
#'
#' @examples
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
#' # Create a list of distance matrices, one for each stratum
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#'
#' @import stats
#' @import MASS
dist_mahal<-function(X, z, st, rank_based = FALSE){
  stopifnot(all((1:length(unique(st)))==sort(unique(st))))
  if (!is.factor(st)) {st <- as.factor(st)}
  Xtemp <- X
  if (is.null(rownames(Xtemp))) {
    rownames(Xtemp) <- 1:nrow(Xtemp)
  }
  if (is.data.frame(X)) X<-as.matrix(X)
  rownames(X) <- rownames(Xtemp)
  colnames(X) <- colnames(Xtemp)

  nobs<-dim(X)[1]
  stopifnot(nobs == length(z))
  if (rank_based) {
    for (j in 1:ncol(X)) {
      X[, j] <- rank(X[, j])
    }
  }
  cv<-cov(X)
  if (rank_based) {
    vuntied <- var(1:nobs)
    rat <- sqrt(vuntied/diag(cv))
    cv <- diag(rat) %*% cv %*% diag(rat)
  }

  icv<-MASS::ginv(cv)
  mats <- list()
  for (j in 1:length(levels(st))) {
    nt <- sum(z[st == j])
    nc <- sum(1-z[st == j])
    if (nt == 0 || nc == 0) {
      warning(paste0("Note that there are either no control or no treated individuals in stratum ",
                     j, ".\n Distance matrix is NA for this stratum.\n"))
      mats[[j]] <- NA
    } else {
      idx_t <- rownames(X[st == j & z == 1, , drop = FALSE])
      idx_c <- rownames(X[st == j & z == 0, , drop = FALSE])
      mats[[j]] <- matrix(NA, nt, nc, dimnames = list(idx_t, idx_c))
      Xst <- X[st == j, , drop = FALSE]
      zst <- z[st == j]

      for (i in 1:nt) {
        mats[[j]][i, ]<-mahalanobis(Xst[zst == 0, , drop = FALSE],
                                    Xst[idx_t[i], , drop = FALSE], icv, inverted=TRUE)
      }
    }
  }
  return(mats)
}
