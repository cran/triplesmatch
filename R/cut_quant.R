#' Create strata based on quantiles of a score
#'
#' @param v Vector of scores (typically propensity scores)
#' @param q Vector of desired quantiles (between 0 and 1) at which to cut the strata
#' @param int Boolean whether to return strata as integers. Default is `TRUE`
#'
#' @return Vector of strata
#' @export
#'
#' @examples
#' cut_quant(1:9, c(1/3, 2/3), int = TRUE)

cut_quant<-function(v,q,int=TRUE){
  qu<-quantile(v,c(0,q,1))
  if (int) {
    as.integer(cut(v,qu,include.lowest=TRUE))
  } else {
    cut(v,qu,include.lowest=TRUE)
  }
}
