#' Formats the triples match for input to other functions
#'
#' @param m The triples match, typically the output of [triples()]
#' @param y A vector or matrix containing the outcome(s). If there are multiple
#' outcomes, there should be a column for each. The rows (or vector elements) correspond to units.
#' The vector should have names or the matrix should have rownames specifying the unit name.
#' These unit names should correspond to those used in `m`
#' @param type Either `wide` or `long`. `wide` formats the match for input to [sentrip()] whereas
#' `long` formats the match for input to [infsentrip()].
#' `wide` can only be used if there is
#' exactly one outcome in `y`. `wide` creates a matrix `ymat` of outcomes with a row
#' corresponding to each triple and three columns corresponding to the units in the triple.
#' It also creates a logical vector `treated1` stating whether the first unit in
#' the corresponding row of `ymat` is the one
#' treated in the triple or not (in which case it would be the one control in the triple).
#' `long` creates a list of three elements: `y`, `z`, and `mset`. Each of these elements
#' is a vector with one element corresponding to each unit in the triples match. `y` is the
#' outcome, `z` is 1 if treated and 0 otherwise, and `mset` is the number of the triple
#' to which this unit belongs
#'
#' @return A list containing either two elements `ymat` and `treated1` if `type == "wide"`
#' or three elements `y`, `z`, and `mset` if `type == "long"`
#' @export
#'
#' @examples
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
#' # Wide version only works with a single outcome
#' ywide <- formattrip(m = triplesm, y = y[, 1], type = "wide")

formattrip <- function(m, y, type = "wide") {
  if (is.list(m) & names(m)[1] == "m") {
    m <- m$m
  }
  if(!all(c("triple", "itreated", "jcontrol", "kthird", "nOfTreated") %in% names(m))) {
    stop('`m` must contain columns for "triple", "itreated", "jcontrol", "kthird", and "nOfTreated"')
  }
  stopifnot(type %in% c("wide", "long"))
  stopifnot(is.vector(y) | is.matrix(y))
  if ((is.vector(y) & is.null(names(y))) | is.matrix(y) & is.null(rownames(y))) {
    stop("`y` must be a vector or matrix with names or rownames corresponding to those
         used in the `itreated`, `jcontrol`, and `kthird` columns of `m`.\n")
  }

  if (type == "wide") {
    if (!is.vector(y)) {
      stop("`y` must be a vector if `type == wide`.\n")
    }
  ymat <- matrix(NA, nrow = nrow(m), ncol = 3)
  treated1 <- rep(NA, nrow(m))
  for (m_no in 1:nrow(m)) {
    if (m[m_no, "nOfTreated"] == 1) {
      ymat[m_no, ] <- y[unlist(m[m_no, c("itreated", "jcontrol", "kthird")])]
      treated1[m_no] <- TRUE
    } else {
      ymat[m_no, ] <- y[unlist(m[m_no, c("jcontrol", "itreated", "kthird")])]
      treated1[m_no] <- FALSE
    }
  }
  return(list(ymat = ymat, treated1 = treated1))
  }

  if (type == "long") {
    if (is.vector(y)) {
      y <- matrix(y, ncol = 1, dimnames = list(names(y), NULL))
    }
    matched <- rownames(y) %in% c(m$itreated, m$jcontrol, m$kthird)

    trip <- data.frame("unit" = rep(NA, 3 * nrow(m)), "triple" = rep(NA, 3 * nrow(m)),
                       "z" = rep(1, 3 * nrow(m)))
    trip[1:nrow(m), 1:2] <- m[, c("itreated", "triple")]
    trip[nrow(m) + 1:nrow(m), 1:2] <- m[, c("jcontrol", "triple")]
    trip[nrow(m) + 1:nrow(m), 3] <- 0
    trip[2*nrow(m) + 1:nrow(m), 1:2] <- m[, c("kthird", "triple")]
    trip[2*nrow(m) + 1:nrow(m), 3][m$nOfTreated == 1] <- 0
    rownames(trip) <- trip$unit
    trip <- trip[rownames(y), ]

    return(list(y = y[matched, ], z = trip$z[matched], mset = trip$triple[matched]))

  }
}
