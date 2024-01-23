#' assess_duplicates
#'
#' Function to assess and resolve elements with multiple higher
#' classifications in a tgraph object. Assessment is performed
#' based on the topology of the graph they form. Linear paths
#' (i.e. two totally separate paths diverging from the a shared
#' node), rings (divergent paths which only reunite at the
#' highest rank in the tgraph) or more than two divergent paths
#' are treated as distinct. If not any of these cases, the
#' distance between the focal element and the reunion of the
#' divergent paths, along with their subtopologies are assessed
#' and a consensus or preferred path based on the frequency of
#' each path in the tgraph or their completeness returned, or
#' the element judged as having multiple distinct
#' classifications
#' @param x A tgraph object
#' @param node character vector of elements with multiple higher
#' classifications in x, or a tvertseq object with those same
#' elements as focal
#' @param mode The rule to be used in selecting between multiple
#' higher classifications. It is possible for the most complete
#' pathway to also be the most frequent
#' @param jump The maximum number of levels between the point of
#' divergence and the point of reunion (if present) for a given
#' path, below which the divergence will be taken as conflicting
#' @param plot A logical speciying if the divergent paths should
#' be plotted
#' @return A list with as many items as elements with multiple
#' classifications, each recording the assessment for a given
#' element
#' @import igraph
#' @import graphics

assess_duplicates <- function(x, node, mode = c("frequency", "completeness"), jump = 3,
                              plot = FALSE) {

  if(!exists("x")) {
    stop("Please supply a tgraph object")
  }
  if(!(inherits(x, "tgraph"))) {
    stop("Please supply a tgraph object")
  }
  if(!exists("node")) {
    stop("Please supply one or more duplicate taxa for testing")
  }
  if(!is.character(node)) {
    stop("Please supply one or more duplicate taxa for testing")
  }
  node <- which(igraph::V(x$taxa)$name %in% node)

  # set frequency as the default decision, if not overwritten in the arguments
  if(length(mode) == 2) {
    mode <- "frequency"
  }
  n <- max(x$ranks)

  # internally define LCA function as is very small
  lca <- function(graph, ...) {
    dots = c(...)
    path = igraph::ego(graph, order = length(igraph::V(graph)), nodes = dots, mode = "out")
    max(Reduce(intersect, path))
  }

  # classify topologies and return list
  ol <- list()
  for(i in 1:length(node)) {

    nd <- node[i]
    # set up return object (out[[2]] is initially set to the failsafe option)
    out <- list()
    out[[1]] <- NA
    out[[2]] <- "Unclassified topology (check that higher level divergences have been resolved first)"
    out[[3]] <- NA
    out[[4]] <- NA
    names(out) <- c("tax", "dec", "prop_del", "prop_add")

    # get complete subgraph for the taxon in question
    paths <- igraph::ego(graph = x$taxa, order = n, nodes = nd, mode = "out")[[1]]
    test <- igraph::induced_subgraph(graph = x$taxa, vids = paths)
    out[[1]] <- paths[1]
    keep_going <- TRUE

    # ranks of the 2-degree nodes ([-1] removes the root vertex)
    vrank <- igraph::V(test)[igraph::ego(test, order = 1, nodes = paths[1]$name)[[1]]]$rank[-1]
    # if there are more than two divergent nodes, the taxon is treated as highly unstable and skipped
    if(length(vrank) > 2) {
      out[[2]] <- paste("Classification indistinct (3+ alternative classifications)")
      keep_going <- FALSE
    }

    # if the subgraph forms is a straight, segmented line (the line test), this means two totally separate classifications,
    # so likely different and so is skipped (topological requirements are n(V) = n(E) + 1, 2 vertex degrees = 1, the rest = 2)
    # i.e. sum if all degrees = 2 would equal 2 * n(V), so having two 1-degree vertices makes sum(degree(V)) == (2 * n(V)) - 2
    if(keep_going & length(igraph::V(test)) == (length(igraph::E(test)) + 1) & sum(igraph::degree(test)) == ((2 * length(igraph::V(test))) - 2)) {
      out[[2]] <- paste("Classification retained (linear divergence)")
      keep_going <- FALSE
    }

    if(keep_going) {
      # as the graph is not a segmented line, there are only two divergence pathways which rejoin at the some point, so get the number,
      # name and rank of the next rank they share
      int <- lca(graph = test, paths[[2]]$name, paths[3]$name)
      v_int <- igraph::V(test)[int]
      n_int <- igraph::V(test)[int]$name
      r_int <- igraph::V(test)[int]$rank
      # if one of the pathways is a direct jump (the jump test) to the shared node with the other pathway, lower classifications have likely been accidentally missed
      # so the jump path is cut, leaving the longer more complete path. This case could theoretically include circles [i.e. a genus to kingdom jump] so it are dealt
      # with before the circle case is tested. (topological requirement means that the jump pathway length is 0)
      p1 <- igraph::all_simple_paths(test, from = paths[2]$name, to = v_int$name)
      p2 <- igraph::all_simple_paths(test, from = paths[3]$name, to = v_int$name)
      if(any(length(p1) == 0 | length(p2) == 0)) {
        out[[2]] <- paste("Reclassification proposed (ranks skipped)")
        out[[3]] <- igraph::E(x$taxa)[paths[1]$name %--% n_int] # record deletions
        keep_going <- FALSE
      }

      # if the subgraph is a full circle (the circle test), the segmented pathways meet only at the top level of the graph (the kingdom level), so
      # they are likely separate and this case is skipped (topological requirements are > 2 vertices, n(V) = n(E), all degree(V) = 2)
      #if(all(igraph::degree(test) == 2) & length(igraph::V(test) == length(igraph::E(test)) & length(igraph::V(test)) > 2)) {
      #  out[[2]] <- paste("Classification retained (circular divergence)")
      #  keep_going <- FALSE
      #}

      # if there are only 4 steps in 2-order neighbourhood of the focal node, the two divergent nodes immediately share a common parent
      # (the immediate parent test) and so are either equivalent (equal rank), or separate parts of the same pathway (different rank)
      if(keep_going) {

        paths2 <- igraph::ego(x$taxa, order = 2, nodes = paths[1], mode = "out")
        if(length(paths2[[1]]) == 4) {

          # if the nodes are of equal rank, delete the lower weighted edge
          if(vrank[1] == vrank[2]) {
            vl <- (igraph::E(x$taxa)[igraph::incident(x$taxa, v = paths[1], mode = "out")])[which.min(igraph::E(x$taxa)[igraph::incident(x$taxa, v = paths[1], mode = "out")]$freq)]
            # delete that edge from the graph, leaving the final classification
            out[[2]] <- paste("Reclassification proposed (equivalent paths)")
            out[[3]] <- vl # record deletions
            # otherwise, combine the different levels into a single pathway
          } else {
            # delete the link between the lower ranked vertex and the parent, and between the higher ranked vertex and the root
            vl <- c((paths[2:3])[which.max(vrank)]$name, n_int)
            vh <- c((paths[2:3])[which.min(vrank)]$name, paths[1]$name)
            out[[2]] <- paste("Reclassification proposed (combined paths)")
            out[[3]] <- igraph::E(x$taxa)[vl[1] %--% vl[2]]
            out[[3]] <- igraph::E(x$taxa)[c(vl[1], vh[2]) %->% c(vl[2], vh[1])] # record deletions
            # link the lower and higher-ranked vertices
            out[[4]] <- igraph::V(x$taxa)[c(paths[2]$name, paths[3]$name)] # record added paths
          }
        } else {

          # if there are only two divergent nodes and the jump, line, circle or immediate parent rules do not apply, the length of the pathways
          # and their frequency are used to decide classification, within a certain rank and jump tolerance
          # get the rank of the focal node
          r_min <- max(V(test)$rank)
          # if the divergent nodes differ by max 1 rank and the next shared node is max 3 levels above the base
          if(diff(vrank) < 2 & (r_min - int) < (jump + 1)) {

            # get the pathway frequencies and lengths
            f1 <- igraph::E(x$taxa)[paths[1] %--% paths[2]]$freq
            f2 <- igraph::E(x$taxa)[paths[1] %--% paths[3]]$freq
            p1 <- length(igraph::all_simple_paths(test, from = paths[2]$name, to = v_int))
            p2 <- length(igraph::all_simple_paths(test, from = paths[3]$name, to = v_int))

            if(mode == "frequency") {
              # deal with equal scores
              if(f1 == f2) {
                if(p1 > p2) {
                  f1 <- f1 + 1
                } else {
                  f2 <- f1 + 1
                }
              }
              vl <- c(paths[2], paths[3])[which.min(c(f1, f2))]
              vl <- c(paths[1], vl)
              out[[2]] <- paste("Reclassification proposed (frequency)")
              out[[3]] <- igraph::E(x$taxa)[vl[1] %--% vl[2]] # record deletions
            }
            if(mode == "completeness") {
              # deal with equal scores
              if(p1 == p2) {
                if(f1 > f2) {
                  p1 <- p1 + 1
                } else {
                  p2 <- p1 + 1
                }
              }
              vl <- c(paths[2], paths[3])[which.min(c(p1, p2))]
              vl <- c(paths[1], vl)
              out[[2]] <- paste("Reclassification proposed (completeness)")
              out[[3]] <- igraph::E(x$taxa)[vl[1] %--% vl[2]] # record deletions
            }

            # otherwise keep the alternative classifications
          } else {
            out[[2]] <- paste("Classification retained (partial divergence)")
          }
        }
      }
    }
    ol[[i]] <- out
  }
  flag <- unlist(lapply(ol, function(y) y$dec))
  if("Unclassified topology (check that higher level divergences have been resolved first)" %in% flag) {
    warning("Unclassified topology encountered")
  }
  return(ol)
}
