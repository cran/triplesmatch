#' Create a triples match for a single stratum
#'
#' @param cost Matrix of distances between treated (rows) and control (columns) units within the stratum
#' @param z Vector of treatment assignments for units within the stratum (0 for control, 1 for treated)
#' @param solver Solver to use for the network problem. Either 'rrelaxiv' or 'rlemon'.
#' 'rrelaxiv' can be downloaded from "https://github.com/josherrickson/rrelaxiv/"
#'
#' @return Named list with three elements: `m` contains the triples match. This is in the form of a
#'  data.frame with number of rows equal to the number of triples and 7 columns specifying the
#'  match number, the names of the three units within the match, the costs of the two
#'  treated-control pairs within the match, and the number of treated units. `obj`
#'  contains the total objective from the network optimization and `bound` contains
#'  a loose lower bound on the objective of the optimal match.
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' n <- 20
#' x <- rnorm(n, 0, 1)
#' nt <- floor(n * 0.62)
#' nc <- n - nt
#' z <- c(rep(1, nt), rep(0, nc))
#' dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
#' triples_st(cost = dist, z = z, solver = "rlemon")
#'
#' @import rlemon
#' @import rcbalance
#' @importFrom optmatch fullmatch

triples_st <- function (cost, z, solver = "rrelaxiv")
{
  if (solver == "rrelaxiv" && !requireNamespace("rrelaxiv", quietly = TRUE)) {
    warning("Package 'rrelaxiv' needed if \"solver\" parameter set to \"rrelaxiv\". \n
            It can be downloaded from \"https://github.com/josherrickson/rrelaxiv/\". \n
            Switching solver to 'rlemon' instead. \n")
    solver <- "rlemon"
  }
  # Check inputs
  stopifnot(is.vector(z) & all((z==0)|(z==1)))
  ntreated <- sum(z)
  ncontrol <- sum(1-z)
  if (ntreated==0){
    warning(paste0("A stratum with ", ncontrol, " controls has zero treated individuals.\n",
                   "No matches produced for this stratum.\n"))
    return(NULL)
  }
  if (ncontrol==0){
    warning(paste0("A stratum with ", ntreated, " treated units has zero control individuals.\n",
                   "No matches produced for this stratum.\n"))
    return(NULL)
  }
  stopifnot((dim(cost)[1])==ntreated) # number of rows must match number of treated units
  stopifnot((dim(cost)[2])==ncontrol) # number of columns must match number of control units

  # Save names of units before reordering by treatment status
  if (!is.null(names(z))) {
    old_names <- names(z)
  } else if (!is.null(rownames(cost)) & !is.null(colnames(cost))) {
    old_names <- rep(NA, length(z))
    old_names[z==1] <- rownames(cost)
    old_names[z==0] <- colnames(cost)
  } else {
    old_names <- 1:length(z)
  }
  # Reorder cost matrix if not in same order as z
  cost <- cost[old_names[z == 1], ]
  cost <- cost[, old_names[z == 0]]

  # Put treated first, then number 1,...,Sample Size
  o <- order(1-z)
  z <- z[o]
  names(z) <- 1:length(z)
  rownames(cost) <- 1:nrow(cost)
  colnames(cost) <- (nrow(cost) + 1):length(z)

  # If all triples are of one type, we will use optmatch except in special circumstances
  use_optmatch <- TRUE

  # Make sure the number of treated units is greater than 0 and
  #     less than or equal to twice the number of controls
  mtk <- pmin(ntreated, 2*ncontrol) # Number of treated to be used
  if (mtk < ntreated) {
    nt_drop <- ntreated - mtk
    warning(paste0("In one stratum, treated individuals exceed twice the controls. \n",
                   nt_drop, " treated will not be matched. \n"))
  } else if (mtk %% 2 == 1 & (floor(mtk / 2) + 2) > ncontrol  ) {
    # If the number of treated units is almost twice the number of controls and is odd,
    #    the last treated unit will need to be matched to two controls.
    #    If this is not possible because there is only one control remaining,
    #    we need to drop both one treated and one control unit.
    # Because we need to drop one of each, we also cannot use optmatch
    mtk <- mtk - 1
    nt_drop <- 1
    use_optmatch <- FALSE
    warning(paste0("In one stratum, 1 treated individual will not be matched due to \n
                   being odd and almost exceeding twice the number of controls.\n"))
  } else {
    nt_drop <- 0
  }

  # bk is the number of triples that will have two treated units
  bk <- max(0, ceiling(((2*mtk)-ncontrol)/3))
  # flow is the number of matches that will need to be made
  nmatches <- mtk-bk
  if (nmatches<0) stop("The number of matches to be made, m_tk-b_k, most be nonnegative.")

  # If all matches will be 2T to 1C but 1 control needs to be dropped due to needing a multiple of 3,
  #   we cannot force this with optmatch
  if (bk == nmatches & ncontrol > nmatches) {
    use_optmatch <- FALSE
  }

  m <- NULL
  # If all matches are of one type and no exceptions occurred, then use optmatch
  if (use_optmatch & (bk == 0 | bk == nmatches)) {
    if (bk == 0) {
      # All matches are 1T:2C and some extra controls will be dropped
      cont <- 2
      nOfTreated <- 1
      omit.fraction <- (ncontrol - 2 * ntreated) / ncontrol
    } else {
      # All matches will be 2T to 1C and some extra treated will be dropped
      cont <- 1/2
      nOfTreated <- 2
      omit.fraction <- - nt_drop/ntreated
    }

    if (is.null(m)) {
      # Run match and return results
      m <- optmatch::fullmatch(cost, min.controls = cont, max.controls = cont,
                               omit.fraction = omit.fraction, data = z)
    }

    # Summary table of matched triples
    ot <- data.frame(matrix(NA, nrow = length(levels(m)), ncol = 7,
                            dimnames = list(NULL, c("triple", "itreated", "jcontrol",
                                                    "kthird", "costStep1", "costStep2", "nOfTreated"))))
    for (m_no in 1:length(levels(m))) {
      m_id <- levels(m)[m_no]
      m_units <- names(m[!is.na(m) & m == m_id])
      t_units <- m_units[m_units %in% rownames(cost)]
      if (t_units[1] != "FakeT") {
        c_units <- m_units[m_units %in% colnames(cost)]
        ot[m_no, "triple"] <- m_no
        ot[m_no, c("itreated", "jcontrol")] <- as.numeric(c(t_units[1], c_units[1]))
        ot[m_no, "costStep1"] <- cost[t_units[1], c_units[1]]
        if (nOfTreated == 2) {
          ot[m_no, "kthird"] <- as.numeric(t_units[2])
          ot[m_no, "costStep2"] <- cost[t_units[2], c_units[1]]
        } else {
          ot[m_no, "kthird"] <- as.numeric(c_units[2])
          ot[m_no, "costStep2"] <- cost[t_units[1], c_units[2]]
        }
        ot[m_no, "nOfTreated"] <- nOfTreated
      }
    }

    # Full match is equal to triples match in this case
    full_bound <- sum(ot$costStep1) + sum(ot$costStep2)

  } else {

    # Step 1 network
    net1temp <- step1network(ntreated = ntreated, ncontrol = ncontrol,
                             nmatches = nmatches, cost = cost)
    net1 <- net1temp$net1

    # Run the minimum cost flow algorithm on the network
    res1 <- rcbalance::callrelax(net1, solver=solver)
    if ((res1$crash==1)) stop("First step crashed")
    if ((!res1$feasible)) stop("First step infeasible")

    # Summary table of matched pairs from first match
    tb <- data.frame("itreated" = net1$startn, "jcontrol" = net1$endn,
                     "used" = res1$x, "costStep1" = net1temp$cost01)[1:net1temp$nedges, ]
    tb <- tb[tb$used==1,]
    tb$triple <- 1:nmatches

    # Step 2 network
    net2temp <- step2network(tb = tb, ntreated = ntreated, ncontrol = ncontrol,
                             nt_drop = nt_drop, cost = cost)
    net2 <- net2temp$net2

    # Run the minimum cost flow algorithm on the network
    res2 <- rcbalance::callrelax(net2, solver=solver)
    if ((res2$crash==1)) stop("Second step crashed")
    if ((!res2$feasible)) stop("Second step infeasible")

    # Summary table of matched triples combining the two steps
    tb2 <- data.frame("triple" = net2$startn, "used" = res2$x,
                      "costStep2" = net2temp$cost02)[1:net2temp$nedges, ]
    tb2$kthird <- net2temp$ids
    tb2 <- tb2[tb2$used==1, ]
    tb2 <- tb2[order(tb2$triple), ]
    ot <- cbind(tb, tb2)
    ot <- ot[, c("triple", "itreated", "jcontrol", "kthird", "costStep1", "costStep2")]
    ot$nOfTreated <- 1 + (ot$kthird<=ntreated)

    # Calculate bound on the total cost based on the corresponding full match
    full_bound <- calc_bound(ncontrol = ncontrol, nt_drop = nt_drop, nmatches = nmatches,
                             cost = cost, z = z, ot = ot)

  }

  # Put back in original order
  ot$itreated <- old_names[o[ot$itreated]]
  ot$jcontrol <- old_names[o[ot$jcontrol]]
  ot$kthird <- old_names[o[ot$kthird]]
  un <- length(unique(c(ot$itreated, ot$jcontrol, ot$kthird)))
  if (un!=(3*nmatches)) warning(paste("Issue: 3*nmatches != unique matched individuals"))

  obj <- sum(ot$costStep1) + sum(ot$costStep2)

  return(list(m = ot, obj = obj, bound = full_bound))
}



