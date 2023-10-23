#' Convert outcome to aberrant ranks
#'
#' Replaces non-aberrant responses by 0 and ranks the aberrant responses by severity.
#' The more aberrant responses have the highest ranks.
#'
#' This function serves the same function as [aberrantscores()] but takes inputs in
#' the `long` format instead of the `wide` format (see [formattrip()] for a
#' description of the two formats, their uses, and their creation).
#'
#' This can be useful for creating a column of `sc` in [infsentrip()] if the
#' aberrant rank test is desired for that variable.
#'
#'
#' @inheritParams aberrantscores
#' @param y Vector of outcomes. Length is equal to the number of units
#' @param tau The null treatment effect to be subtracted from all treated units
#' before aberrant ranking commences. If `tau != 0`, then `z` is required
#' @param z Vector with length equal to that of y. Each element specifies whether a unit is treated (1) or not (0)
#'
#' @return Vector of aberrant ranks corresponding to `y`
#' @export
#' @seealso aberrantscores for the same function with inputs in the wide format instead of long
#' @seealso formattrip for formatting of the triples match into wide or long format
#'
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
#' # Aberrant ranks for primary outcome
#' y[, 1] <- aberrantscoreslong(y[, 1], cutoff = 0.5, cutoff_dir = "greater")

aberrantscoreslong <- function(y, cutoff, cutoff_dir = "less", tau = 0, z = NULL) {
  stopifnot(cutoff_dir %in% c("less", "greater"))
  stopifnot(is.vector(y))

  if (!(tau == 0)) {

    if (is.null(z)) {
      stop("If `tau` is nonzero, `z` vector must be supplied.\n")
    }
    stopifnot(length(y) == length(z))
    stopifnot(all(z %in% c(0,1)))
    stopifnot(sum(z) > 0)
    stopifnot(sum(1-z) > 0)

    # Subtract the null effect from treated units
    y[z == 1] <- y[z==1] - tau
  }

  if (cutoff_dir == "less") {
    aberrant <- y < cutoff
    aberrant_ranks <- rank(-y[aberrant], ties.method = "average")
  } else {
    aberrant <- y > cutoff
    aberrant_ranks <- rank(y[aberrant], ties.method = "average")
  }

  # Replace aberrant outcomes with their ranks
  y[aberrant] <- aberrant_ranks
  # Replace non-aberrant outcomes with 0
  y[!aberrant] <- 0

  return(y)
}
