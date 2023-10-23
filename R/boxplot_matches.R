#' Series of boxplots for a given variable characterizing the triples match
#'
#' @param yname y axis label
#' @param m `m` element of the list returned from `triples()` function containing information
#' about matched individuals
#' @param y Named vector containing variable to plot on the y axis.
#' Names must correspond to the units specified in `m`
#' @param z Vector of treatment indicators. Must be in same order as `y`
#' @return Display containing three sets of boxplots for the propensity score.
#'     First is for all treated vs control units.
#'     Second is for the triples that have one treated unit and two controls.
#'     Third is for the triples that have two treated units and one control.
#' @export
#'
#' @import graphics
#'
#' @examples
#' # Generate some data
#' set.seed(8)
#' n <- 200
#' nt <- floor(n * 0.5)
#' nc <- n - nt
#' x <- c(rnorm(nt, 0, 1), rnorm(nc, 0.6, 1))
#' z <- c(rep(1, nt), rep(0, nc))
#' # Create some strata
#' ps <- glm(z ~ x, family = binomial)$fitted.values
#' ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
#' # Create a distance matrix
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#' # Construct the triples match
#' triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
#' boxplot_matches(m = triplesm$m, y = ps, z = z, yname = "Propensity score")

boxplot_matches<-function(m, y, z, yname = NULL){
  stopifnot(!is.null(names(y)))
  if (!is.factor(z)) {
    z <- factor(z, levels = c(1, 0), labels = c("T", "C"))
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, 3))
  par(mai = c(0.75, 0.5, 0.5, 0.25), oma = c(0, 1.5, 2, 0))
  v1 <- y[m$itreated]
  v2 <- y[m$jcontrol]
  v3 <- y[m$kthird]
  k <- (m$nOfTreated==1)
  mn <- min(y)
  mx <- max(y)
  boxplot(y ~ z, names = c("T", "C"), ylab = yname, cex.main = .9, cex.lab = .9,
          cex.axis = .9, las = 1, main = paste0(length(z), " Unmatched"),
          xlab = "Treatment Group", ylim = c(mn, mx))
  boxplot(v1[k], v2[k], v3[k], names = c("T", "C", "C"), ylab = NULL,
          ylim = c(mn, mx), xlab = "(i, j, k)", cex.main = .9, cex.lab = .9,
          cex.axis = .9, las = 1, main = paste("1 Treated, ", sum(k), "triples"))
  boxplot(v1[!k], v2[!k], v3[!k], names = c("T", "C", "T"), ylab = NULL,
          ylim = c(mn, mx), xlab = "(i, j, k)", cex.main = .9, cex.lab = .9,
          cex.axis = .9, las = 1, main = paste("2 Treated, ", sum(!k), "triples"))
}
