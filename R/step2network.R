#' Second network of the triples match algorithm
#'
#' Creates the network that is used in the second step of the triples match.
#' This network has edges between previously constructed pair matches
#' to the remaining available units.
#'
#' @param tb Data.frame containing a summary from first step.
#' Columns include `itreated`, `jcontrol`, `used`, `costStep1`, `triple`
#' @param ntreated Number of treated units
#' @param ncontrol Number of control units
#' @param nt_drop Number of treated units to be dropped
#' @param cost The cost matrix to be used in the optimization. Number of rows
#' should equal `ntreated` and number of columns should equal `ncontrol`
#'
#' @return Named list with four elements: `net2` contains the network,
#' `cost02` contains the true costs of the edges in the network before
#' they are manipulated for the purposes of the network optimization,
#' `nedges` contains the number of edges between pairs and remaining units,
#' and `ids` contains the names of the remaining individual involved in each edge.
#' The network `net2` is itself a named list, containing five elements:
#' `startn`, `endn`, `ucap`, `cost`, `b`. These describe the starting node,
#' ending node, capacity, cost, and supply of each of the edges in the network.
#' @noRd

step2network <- function(tb, ntreated, ncontrol, nt_drop, cost) {
  stopifnot(nrow(cost) == ntreated)
  stopifnot(ncol(cost) == ncontrol)

  # Network 2 has used pairs vs unused individuals

  nmatches <- nrow(tb)
  ij <- 1:nmatches

  # Create node numbers for unused individuals
  # k contains id numbers for individuals not matched in first match
  # kk assigns them node numbers
  k <- (1:(ntreated+ncontrol))
  k <- k[!is.element(k, c(tb$itreated, tb$jcontrol))]
  kk <- nmatches+1:length(k)

  # Nodes
  # 1 to nmatches are the used pairs, the next length(k) are the unused units (control or treated)
  nodes2 <- c(ij,kk)

  # Add a source and a sink
  source <- length(nodes2)+1
  sinkT <- length(nodes2)+2
  sinkC <- length(nodes2)+3
  nodes2 <- c(nodes2,source,sinkT,sinkC)

  # Flow supply/demand
  b2 <- c(rep(0,length(nodes2)-3),
          nmatches,
          -(ntreated-nt_drop-nmatches), # Treated sink demands the rest of the treated to be used
          -((2*nmatches)-ntreated+nt_drop)) # Control sink demands the rest of the supply

  # Edges
  # Pair-unused edges
  # Start at each pair going to unused individual 1, then repeat for remaining unused:
  #     Pair1 to unused1, pair2 to unused1, ...
  startn2 <- rep(ij,length(kk))
  endn2 <- rep(kk, each = nmatches)
  # Create the costs
  cost2 <- matrix(NA,nmatches,length(kk))
  rownames(cost2) <- ij
  colnames(cost2) <- kk
  for (p in ij){
    whoi <- tb$itreated[p] # Treated in the pair
    whoj <- tb$jcontrol[p] # Control in the pair
    for (u in (kk-nmatches)){
      whok <- k[u] # Unused person to add
      # Cost is always treated-control
      #    Between unused person and whoever has opposite treatment in the pair
      if (whok<=ntreated) cost2[p,u] <- cost[whok, whoj - ntreated]
      else cost2[p, u] <- cost[whoi, whok - ntreated]
    }
  }
  cost2 <- as.vector(cost2)
  # Drop any edges that have infinite cost
  startn2 <- startn2[cost2 < Inf]
  endn2 <- endn2[cost2 < Inf]
  k2 <- rep(k, each = nmatches)
  k2 <- k2[cost2 < Inf]
  cost2 <- cost2[cost2 < Inf]
  n_tc_edges2 <- length(startn2)

  # Source-pair edges (cost 0)
  startn2 <- c(startn2, rep(source, nmatches))
  endn2 <- c(endn2, ij)
  cost2 <- c(cost2, rep(0, nmatches))

  # Unused person - corresponding sink edges (cost 0)
  startn2 <- c(startn2, kk)
  for (u in (kk-nmatches)) {
    whok <- k[u]
    if (whok<=ntreated) endn2 <- c(endn2, sinkT)
    else endn2 <- c(endn2, sinkC)
  }
  cost2 <- c(cost2, rep(0, length(kk)))

  # Add stability increment and convert to integers
  cost02 <- cost2
  cost2 <- (cost2 + 1)/1e-05
  cost2 <- round(cost2)

  # Every edge has capacity 1
  ucap2  <-  c(rep(1, length(startn2)))

  # If none of the edges end at sinkT, remove it and change the number of sinkC accordingly
  if (!sinkT %in% endn2) {
    endn2[endn2 == sinkC] <- sinkC - 1
    nodes2 <- nodes2[-length(nodes2)]
    b2 <- b2[-sinkT]
    sinkC <- sinkC - 1
    sinkT <- NA
  }

  # If none of the edges end at sinkC, remove it
  if (!sinkC %in% endn2) {
    nodes2 <- nodes2[-length(nodes2)]
    b2 <- b2[-sinkC]
    sinkC <- NA
  }
  net2 = list(startn = startn2, endn = endn2, ucap = ucap2, cost = cost2, b = b2)

  return(list(net2 = net2, cost02 = cost02, nedges = n_tc_edges2, ids = k2))
}
