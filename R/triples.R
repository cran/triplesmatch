#' Create a triples match
#'
#' @param cost List of matrices of distances between treated (rows) and control (columns) units within a stratum
#'        with one entry in the list per stratum
#' @param z Vector of treatment assignment (0 for control, 1 for treated)
#' @param st Vector of stratum assignments
#' @param solver Solver to use for the network problem. Either 'rrelaxiv' or 'rlemon'.
#' 'rrelaxiv' can be downloaded from "https://github.com/josherrickson/rrelaxiv/"
#'
#' @return Named list with three elements: `m` contains the triples match. This is in the form of a
#'  data.frame with number of rows equal to the number of triples and 8 columns specifying the
#'  match number, the names of the three units within the match, the costs of the two
#'  treated-control pairs within the match, the number of treated units, and the stratum. `obj`
#'  contains the total objective from the network optimization and `bound` contains
#'  a loose lower bound on the objective of the optimal match.
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
#' # Create a distance matrix
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#' # Construct the triples match
#' triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
#'
#' @import rlemon
#' @import rcbalance
triples <- function (cost, z, st, solver = "rrelaxiv")
{
  if (solver == "rrelaxiv" && !requireNamespace("rrelaxiv",
                                                quietly = TRUE)) {
    warning("Package 'rrelaxiv' needed if \"solver\" parameter set to \"rrelaxiv\". \n It can be downloaded from \"https://github.com/josherrickson/rrelaxiv/\". \n Switching solver to 'rlemon' instead.")
    solver <- "rlemon"
  }
  # Check input
  stopifnot(is.list(cost))
  K <- length(unique(st))
  stopifnot(length(cost) == K)
  stopifnot(length(st)==length(z))
  old_st_names <- sort(unique(st))
  if (!is.factor(st) | !all(old_st_names == 1:K)) {
    st <- factor(st, levels = old_st_names, labels = 1:K)
  }
  for (i in 1:K) {
    stopifnot((dim(cost[[i]])[1])==sum(z[st == i])) # number of rows must match number of treated units
    stopifnot((dim(cost[[i]])[2])==sum(1-z[st == i])) # number of columns must match number of control units
  }
  stopifnot(is.vector(z) & all((z==0)|(z==1)))

  for (i in 1:K) {
    z_ist <- z[st == i]
    if (is.null(names(z_ist)) & !is.null(rownames(cost[[i]])) & !is.null(colnames(cost[[i]]))) {
      names(z_ist)[z_ist==1] <- rownames(cost[[i]])
      names(z_ist)[z_ist==0] <- colnames(cost[[i]])
    }
  }

  if(is.null(names(z))) {
    names(z) <- 1:length(z)
  }

  last_triple <- 0
  m <- NULL

  bound <- 0
  for (i in 1:K) {

    z_ist <- z[st == i]

    if (sum(z_ist) > 0 & length(z_ist) >= 3) {
      if (is.null(rownames(cost[[i]])) | !is.null(colnames(cost[[i]]))) {
        rownames(cost[[i]]) <- names(z_ist)[z_ist == 1]
        colnames(cost[[i]]) <- names(z_ist)[z_ist == 0]
      }

      m_st_temp <- triples_st(cost = cost[[i]], z = z_ist, solver = solver)
      if (!is.null(m_st_temp)) {
        m_st <- m_st_temp$m
        bound <- bound + m_st_temp$bound
        ind_st <- which(st == i)
        m_st$triple <- last_triple + m_st$triple
        m_st$st <- levels(st)[i]
        m <- rbind(m, m_st)
        last_triple <- max(m$triple)
      }
    }
  }

  obj <- sum(m$costStep1) + sum(m$costStep2)

  # Return to original stratum labels
  m$st <- old_st_names[as.numeric(m$st)]

  return(list(m = m, obj = obj, bound = bound))
}

