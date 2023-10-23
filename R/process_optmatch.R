#' Process the output of the optmatch for reporting in simulation studies
#'
#' Of limited interest to most users, this computes information
#' about the quality of a pair or full match and its power. This is used in the simulations in
#' the manuscript "A New Design for Observational Studies Applied to the Study
#' of the Effects of High School Football on Cognition Late in Life" by Brumberg et al.
#' To evaluate simulated power, sensitivity analyses are run using [sensitivityfull::senfm()].
#' This function uses Huber's M-statistic in the hypothesis tests.
#'
#' @inheritParams process_triples_match
#' @param m The resulting match from [optmatch::optmatch()]
#' @param cost_unpen The unpenalized cost matrix
#' @param cost_pen The cost matrix with an infinite penalty added to the distance between
#' any pair of units that are in different strata
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
#'
#' @noRd
process_optmatch <- function(m, cost_unpen, cost_pen, gamma_list, Y, trim = 3) {
  n_matches <- length(levels(m))

  if(is.null(rownames(cost_unpen)) | is.null(colnames(cost_unpen)) |
     is.null(rownames(cost_pen)) | is.null(colnames(cost_pen))) {
    stop("All cost matrices must have row and column names corresponding with the names of `Y` and the names used in the match")
  }
  if(is.null(names(Y))) {
    stop("`Y` must have names matching the rownames and colnames of the cost matrices")
  }

  cost_unpen_sum <- 0 # Sum of unpenalized costs
  cost_pen_all <- NULL # vector of penalized costs, including those that cross strata
  n_crosses <- 0
  n21 <- 0
  n11 <- 0
  n12 <- 0
  ymat <- matrix(NA, nrow = n_matches, ncol = 3)
  treated1 <- rep(NA, n_matches)

  for (m_no in 1:n_matches) {
    m_id <- levels(m)[m_no]
    m_units <- names(m[!is.na(m) & m == m_id])
    t_units <- m_units[m_units %in% rownames(cost_unpen)]
    c_units <- m_units[m_units %in% colnames(cost_unpen)]

    if (length(t_units) == 1) {
      if (length(c_units) == 1) {
        n11 <- n11 + 1
        ymat[m_no, ] <- Y[(c(t_units, c_units, NA))]
      } else {
        n12 <- n12 + 1
        ymat[m_no, ] <- Y[(c(t_units, c_units))]
      }
      treated1[m_no] <- TRUE
    } else {
      n21 <- n21 + 1
      ymat[m_no, ] <- Y[(c(c_units, t_units))]
      treated1[m_no] <- FALSE
    }

    for (t_unit in t_units) {
      for (c_unit in c_units) {
        cost_unpen_sum <- cost_unpen_sum + cost_unpen[t_unit, c_unit]
        if (cost_pen[t_unit, c_unit] > 1000) {
          n_crosses <- n_crosses + 1
        }
        cost_pen_all <- c(cost_pen_all, cost_pen[t_unit, c_unit])

      }
    }
  }

  cost_quantiles <- quantile(cost_pen_all, c( 0.5, 0.75, 0.90, 0.95), na.rm = TRUE)
  names(cost_quantiles) <- c("cost_med", "cost_q75", "cost_q90", "cost_q95")

  quality <- c("n_crosses" = n_crosses,
               "n_t_used" = sum(!is.na(m[names(m) %in% rownames(cost_unpen)])),
               "n_c_used" = sum(!is.na(m[names(m) %in% colnames(cost_unpen)])),
               "n12" = n12,
               "n11" = n11,
               "n21" = n21,
               "cost_subset" = NA,
               "cost_total_unpen" =  cost_unpen_sum,
               "cost_avg_unpen" =  cost_unpen_sum / length(cost_pen_all),
               "cost_avg_pen" =  mean(cost_pen_all),
               cost_quantiles)

  # Power
  reject <- rep(0, length(gamma_list))
  names(reject) <- as.character(gamma_list)
  for (gamma in gamma_list) {
    if (sensitivityfull::senfm(y = ymat, treated1 = treated1, gamma = gamma, inner = 0, trim = trim)$pval <= 0.05) {
      reject[as.character(gamma)] <- 1
    }
  }

  return(list(quality = quality,
              reject = reject))

}
