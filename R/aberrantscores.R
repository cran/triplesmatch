#' Convert outcome to aberrant ranks
#'
#' Replaces non-aberrant responses by 0 and ranks the aberrant responses by severity.
#' The more aberrant responses have the highest ranks.
#'
#' This can be useful for creating `scores` in [sentrip()] for an aberrant rank test.
#'
#' @inheritParams sentrip
#' @param ymat A matrix of outcomes. Rows correspond to matched triples and
#' the three columns correspond to the three units in the match. The first unit
#' is the one treated unit if `treated1 == TRUE` or the one control unit if
#' `treated1 == FALSE`. The other two columns contain the remaining two units in
#' the match. These are control units if `treated1 == TRUE` or treated units if
#' `treated1 == FALSE`. This can easily be created from the triples match using
#' the [formattrip()] function with `type == "wide"`
#' @param cutoff The cutoff for whether an outcome is aberrant. Any outcome more extreme
#' then this cutoff will be considered aberrant
#' @param cutoff_dir Either `less` or `greater`, indicating whether outcomes should be
#' aberrant if they are less than the `cutoff` or greater than the `cutoff`, respectively
#' @param tau The null treatment effect to be subtracted from all treated units
#' before aberrant ranking commences. If `tau != 0`, then `treated1` is required
#'
#' @return A matrix similar to `ymat` in all regards other than the outcomes being
#' converted to aberrant ranks
#' @export
#' @seealso aberrantscoreslong for the same function with inputs given in the long format as opposed to the wide format
#' @seealso formattrip for formatting the triples match into long or wide format
#'
#' @examples
#' # Generate some data
#' set.seed(246)
#' n <- 30
#' x <- rnorm(n, 0, 3)
#' nt <- floor(n * 0.5)
#' nc <- n - nt
#' z <- c(rep(1, nt), rep(0, nc))
#' # Create a distance matrix, everything in one stratum
#' dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
#' # Create the triples match
#' triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")
#' # Create an outcome
#' y <- 1:40
#' # Give the outcome some random unit names
#' names(y) <- sample(1:40)
#' # Reformat the triples match
#' ywide <- formattrip(m = triplesm, y = y, type = "wide")
#' # Turn the outcome into scores, in this case aberrant ranks
#' ab <- aberrantscores(ywide$ymat, 15, cutoff_dir = "less", tau = 0, treated1 = NULL)
#' # Conduct a one-sided hypothesis test with a bias of gamma = 1.25
#' sentrip(scores = ab, treated1 = ywide$treated1, gamma = 1.25, alternative = "greater")
#'

aberrantscores <- function(ymat, cutoff, cutoff_dir = "less", tau = 0, treated1 = NULL) {
  stopifnot(cutoff_dir %in% c("less", "greater"))
  stopifnot(is.matrix(ymat) | is.data.frame(ymat))
  if (is.data.frame(ymat)) {
    ymat <- as.matrix(ymat)
  }
  stopifnot(all(!is.na(as.vector(ymat[, 1:2]))))

  if (!(tau == 0)) {

    if (is.null(treated1)) {
      stop("If `tau` is nonzero, `treated1` vector must be supplied.\n")
    }
    stopifnot(is.logical(treated1))
    stopifnot((dim(ymat)[1]) == length(treated1))

    if (sum(treated1) > 0)   # Are there any matches with 1 treated?
      ymat[treated1, 1] <- ymat[treated1, 1] - tau   # Subtract the null effect from treated units
    if (sum((1 - treated1)) > 0)   # Are there any matches with 1 control?
      ymat[!treated1, -1] <- ymat[!treated1, -1] - tau   # Subtract the null effect from non-control units
  }

  y_nrow <- nrow(ymat)
  y_ncol <- ncol(ymat)
  y <- as.vector(ymat)

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

  # Convert back to matrix
  ymat_out <- matrix(data = y, nrow = y_nrow, ncol = y_ncol, byrow = FALSE,
                     dimnames = dimnames(ymat))

  return(ymat_out)
}
