#' Process the output of the triples match for reporting in simulation studies
#'
#' Of limited interest to most users, this computes information
#' about the quality of a triples match and its power. This is used in the simulations in
#' the manuscript "A New Design for Observational Studies Applied to the Study
#' of the Effects of High School Football on Cognition Late in Life" by Brumberg et al.
#' To evaluate simulated power, sensitivity analyses are run using [sensitivityfull::senfm()].
#' This function uses Huber's M-statistic in the hypothesis tests.
#'
#'
#' @param m The output of [triples()]
#' @param gamma_list The gammas at which to conduct the sensitivity analyses
#' @param Y The outcome to be used for the hypothesis tests
#' @param trim The amount of trimming to use in Huber's M-statistic for the hypothesis tests
#'
#' @return A named list with two elements: `quality` is a named vector describing the
#' number of times strata are crossed in the triples match (this is always 0),
#' the number of treated and control units used in the match, the number of 1T:2C,
#' 1T:1C, and 2T:1C matches constructed (1T:1C is always 0), the cost of the
#' first step of the triples algorithm, the total cost, the average cost across
#' matches, the average penalized cost across matches where the penalty is for
#' crossing strata (this is equal to the unpenalized cost since no stratum
#' crossings occur), and the quantiles of the penalized cost; `reject` is
#' a named vector with a 0 or 1 describing whether the hypothesis test was rejected
#' at each gamma value in `gamma_list`.
#'
#' @noRd
process_triples_match <- function(m, gamma_list, Y, trim = 3) {
  if (is.null(names(Y))) {
    names(Y) <- 1:length(Y)
  }

  triple_cost_all <- c(m$costStep1, m$costStep2)
  cost_quantiles <- quantile(triple_cost_all, c( 0.5, 0.75, 0.90, 0.95), na.rm = TRUE)
  names(cost_quantiles) <- c( "cost_med", "cost_q75", "cost_q90", "cost_q95")

  avg_cost = mean(c(m$costStep1, m$costStep2))
  total_cost = sum(c(m$costStep1, m$costStep2))

  # power
  ymat <- matrix(NA, nrow = nrow(m), ncol = 3)
  treated1 <- rep(NA, nrow(m))

  for (m_idx in 1:nrow(m)) {
    if (m[m_idx, "nOfTreated"] == 1) {
      ymat[m_idx, ] <- Y[unlist(m[m_idx, c("itreated", "jcontrol", "kthird")])]
      treated1[m_idx] <- TRUE
    } else {
      ymat[m_idx, ] <- Y[unlist(m[m_idx, c("jcontrol", "itreated", "kthird")])]
      treated1[m_idx] <- FALSE
    }
  }

  reject <- rep(0, length(gamma_list))
  names(reject) <- as.character(gamma_list)
  for (gamma in gamma_list) {
    if (sensitivityfull::senfm(y = ymat, treated1 = treated1, gamma = gamma, inner = 0, trim = trim)$pval <= 0.05) {
      reject[as.character(gamma)] <- 1
    }
  }

  return(list(quality = c("n_crosses" = 0,
                          "n_t_used" = sum(m$nOfTreated),
                          "n_c_used" = 3 * nrow(m) - sum(m$nOfTreated),
                          "n12" = sum(m$nOfTreated == 1),
                          "n11" = 0,
                          "n21" = sum(m$nOfTreated == 2),
                          "cost_subset" = mean(m$costStep1),
                          "cost_total_unpen" = total_cost,
                          "cost_avg_unpen" = avg_cost,
                          "cost_avg_pen" = avg_cost,
                          cost_quantiles),
              reject = reject
  ))
}
