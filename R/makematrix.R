#' Separate covariate matrix by type of match
#'
#' @param X Covariate matrix
#' @inheritParams boxplot_matches
#'
#' @return data.frame with the first five columns: `triple` (match number),
#' `type` (whether the unit is from `itreated`, `jcontrol`, or `kthird` in the match),
#' `index` (unit name),
#' `z` (treatment indicator),
#' `N` (number of treated units in the match), and
#' additional columns for the covariates supplied in `X`
#' @noRd

makematrix <- function(X, z, m){
  stopifnot(!is.null(rownames(X)))
  stopifnot(!is.null(colnames(X)))
  stopifnot(!is.null(names(z)))
  if (is.vector(X)) X <- matrix(X, length(X), 1)
  triple <- m$triple
  type <- rep(1, dim(m)[1])
  index <- m$itreated
  z1 <- z
  z <- z1[index]
  x <- X[index, ]
  N <- m$nOfTreated
  o <- data.frame(triple, type, index, z, N, x)
  type <- rep(2, dim(m)[1])
  index <- m$jcontrol
  z <- z1[index]
  x <- X[index, ]
  o <- rbind(o, data.frame(triple, type, index, z, N, x))
  type <- rep(3, dim(m)[1])
  index <- m$kthird
  z <- z1[index]
  x <- X[index, ]
  o <- rbind(o, data.frame(triple, type, index, z, N, x))
  o <- o[order(o$triple, o$type), ]
  o
}
