#' First network of the triples match algorithm
#'
#' Creates the network that is used in the first step of the triples match.
#' This network has edges between treated and control units.
#'
#' @param ntreated Number of treated units
#' @param ncontrol Number of control units
#' @param nmatches Number of pairs to be created
#' @param cost Cost matrix. The number of rows should be equal to `ntreated` and
#' the number of columns should be equal to `ncontrol`
#'
#' @return Named list with three elements: `net1` contains the network,
#' `cost01` contains the true costs of the edges in the network before
#' they are manipulated for the purposes of the network optimization and
#' `nedges` contains the number of edges between treated and control units.
#' The network `net1` is itself a named list, containing five elements:
#' `startn`, `endn`, `ucap`, `cost`, `b`. These describe the starting node,
#' ending node, capacity, cost, and supply of each of the edges in the network.
#' @noRd
#'
step1network <- function(ntreated, ncontrol, nmatches, cost) {
  # Nodes
  # 1 to ntreated are treated, ntreated+1 to ntreated+ncontrol are controls
  nodes1 <- 1:(ntreated+ncontrol)

  # Add a source and a sink
  source <- length(nodes1)+1
  nodes1 <- c(nodes1, source)
  sink <- length(nodes1)+1
  nodes1 <- c(nodes1, sink)

  # Flow supply/demand
  # Source has supply of nmatches and sink has demand of nmatches
  b1 <- c(rep(0, ntreated+ncontrol), nmatches, -nmatches)

  # Edges
  # Treated-control edges (cost equal to distance metric between the treated and control)
  #     In order: T1 to C1, T2 to C1, ..., T(I-1) to CJ, TI to CJ
  startn1 <- rep(1:ntreated, ncontrol)
  endn1 <- rep((ntreated+1):(ntreated+ncontrol), each = ntreated)
  cost1 <- as.vector(cost)
  # Drop any edges that have infinite cost
  startn1 <- startn1[cost1 < Inf]
  endn1 <- endn1[cost1 < Inf]
  cost1 <- cost1[cost1 < Inf]
  n_tc_edges <- length(startn1)

  # Source-treated edges (cost 0)
  startn1 <- c(startn1, rep(source, ntreated))
  endn1 <- c(endn1, 1:ntreated)
  cost1 <- c(cost1, rep(0, ntreated))

  # Control-sink edges (cost 0)
  startn1 <- c(startn1, (ntreated+1):(ntreated+ncontrol))
  endn1 <- c(endn1, rep(sink, ncontrol))
  cost1 <- c(cost1, rep(0, ncontrol))

  # Adding 1 to all costs helps with stability (see section 4.1 of Hansen & Klopfer 2006)
  # Cost inputs to callrelax should be integers --
  #    Divide by a common denominator, aka a tiny number, and then round
  cost01 <- cost1 # Save the original costs for reporting later
  cost1 <- (cost1 + 1)/1e-05
  cost1 <- round(cost1)

  # Every edge has capacity 1
  ucap1  <-  c(rep(1, length(startn1)))

  net1 = list(startn = startn1, endn = endn1, ucap = ucap1, cost = cost1, b = b1)
  return(list(net1 = net1, cost01 = cost01, nedges = n_tc_edges))
}
