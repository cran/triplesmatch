#' Boxplots of pairwise differences in triples match
#'
#' Make boxplots of treated - control pair differences before matching, for the two types of triples, and weighted across triples
#'
#' @inheritParams boxplot_matches
#' @return Boxplots with treated minus control pair differences for the specified covariate.
#'     Boxplots are show for before matching, for the matches with 1 treated individual,
#'     for the matches with 2 treated individuals, and for the weighted combination that
#'     duplicates the differences for the matches with two treated individuals
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
#' # Create some strata
#' ps <- glm(z ~ x, family = binomial)$fitted.values
#' ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
#' # Create a distance matrix
#' dist <- dist_mahal(data.frame(x = x), z, ps_st)
#' # Construct the triples match
#' triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
#' boxplot_diffs(m = triplesm$m, y = ps, z = z, yname = "Propensity score")

boxplot_diffs <- function(m, y, z, yname=NULL){
  y1 <- y[m$itreated]
  y2 <- y[m$jcontrol]
  y3 <- y[m$kthird]
  k <- (m$nOfTreated==1)
  mn <- min(y)
  mx <- max(y)
  tc1 <- c(y1[k]-y2[k], y1[k]-y3[k])
  tc2 <- c(y1[!k]-y2[!k], y3[!k]-y2[!k])
  tc <- c(tc2, tc2, tc1)
  tc0 <- as.vector(outer(y[z==1], y[z==0], "-"))

  boxplot(tc0, tc1, tc2, tc, names=c("B", "1", "2", "W"), ylab="", main=yname,
          xlab="", cex.main=.8, cex.lab=.8, cex.axis=.8, las=1)
  abline(h=0)
}
