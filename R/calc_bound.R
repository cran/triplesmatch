#' Calculate lower bound on the triples match objective using full match
#' 
#' This is an auxiliary function, not of interest to users.
#'
#' @param ncontrol Number of controls
#' @param nt_drop Number of treated units to be dropped
#' @param nmatches Number of matches to create
#' @param cost Distance matrix
#' @param z Treatment vector
#' @param ot Summary information generated from triples match
#'
#' @return Numeric bound for the triples match objective
#' @noRd


calc_bound <- function(ncontrol, nt_drop, nmatches, cost, z, ot) {
  
  # In the case of dropping > 1 treated, the triples match = full match, 
  #    so the bound is already taken care of
  
  if (ncontrol - nmatches == 1 & nt_drop == 1) {
    # Need to drop 1 control and 1 treated
    # Optmatch doesn't like this, but we can manipulate it
    # Add two fake controls and two fake treated units
    # The first fake control has cost 0 to all treated units,
    #    and is forced to be matched with the first fake treated
    #    that will have cost 0 to the fake control but infinite to all others.
    #    The full match has the option of dropping one treated unit by adding it to
    #    this pair match with 0 cost. It may choose not to if there is a better solution
    #    that includes all treated units.
    # The second fake control will have cost infinite to all treated
    #    except the second fake treated unit, which will have cost 0.
    #    That fake treated unit will have cost 0 to all real controls and second fake control
    #    and infinite to the first fake control.  The full match has the option of dropping
    #    one control unit by adding it to this pair match with 0 cost. 
    #    It may choose not to if there is a better solution that includes all control units.
    
    # Add column for fake control 1. Cost 0 to all real treated
    cost <- cbind(cost, rep(0, nrow(cost)))
    # Add row for fake treated 1. Cost 0 to fake control 1, infinite for others
    cost <- rbind(cost, c(rep(Inf, ncol(cost) - 1), 0))
    # Add column for fake control 2. Infinite to all real treated and fake treated 1
    cost <- cbind(cost, rep(Inf, nrow(cost)))
    # Add row for fake treated 2. Cost 0 to real controls and fake control 2,
    #     infinite for fake control 1
    cost <- rbind(cost, c(rep(0, ncol(cost) - 2), Inf, 0))
    
    z <- c(z, 0, 1, 0, 1)
    names(z)[length(z) - 3:0] <- c("fakeC1", "fakeT1", "fakeC2", "fakeT2")
    row.names(cost)[nrow(cost) - 1:0] <- c("fakeT1", "fakeT2")
    colnames(cost)[ncol(cost) - 1:0] <- c("fakeC1", "fakeC2")
    
    fullm <- optmatch::fullmatch(cost, min.controls = 1/2, max.controls = 2,
                                 omit.fraction = 0, data = z)
    
  } else  {
    # Try to drop the same amount of controls as in the triples
    # It may drop fewer units, but it will only do this if the objective is lower, 
    #    meaning the bound is more conservative but still valid
    c_drop <- ncontrol - (3 * nmatches - sum(ot$nOfTreated))
    omit_frac <- 0
    if (c_drop > 0) {
      omit_frac <- c_drop / ncontrol
    }
    
    suppressWarnings(
      fullm <- optmatch::fullmatch(x = cost,
                                   min.controls = 1/2, max.controls = 2,
                                   omit.fraction = omit_frac))
  }
  
  if (sum(optmatch::matchfailed(fullm)) > 0) {
    full_bound <- NA
  } else {
    n_matches <- length(levels(fullm))
    
    full_bound <- 0
    for (m_no in 1:n_matches) {
      m_id <- levels(fullm)[m_no]
      m_units <- names(fullm[!is.na(fullm) & fullm == m_id])
      if (all(!c("fakeC1", "fakeT1", "fakeC2", "fakeT2") %in% m_units)) {
        t_units <- m_units[m_units %in% rownames(cost)]
        c_units <- m_units[m_units %in% colnames(cost)]
        for (t_unit in t_units) {
          for (c_unit in c_units) {
            full_bound <- full_bound + cost[t_unit, c_unit]
          }
        }
      }
    }
  }
  
  return(full_bound)
}
