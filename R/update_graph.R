#' update_graph
#'
#' Function to update the structure of a graph, given
#' a set of modification as returned by assess_duplicates
#' @param x a tgraph object to modify
#' @param del A vector of element names or numbers to delete
#' @param add An edge sequence of edges to add to the graph
#' @param changes Alternatively, the output of
#' assess_duplicates, containing proposed deletions and
#' additions
#' @return An updated tgraph object
#' @import igraph

update_graph <- function(x, del = NULL, add = NULL, changes = NULL) {

  if(!exists("x")) {
    stop("Please supply tgraph object for modification")
  }

  if(!is.null(changes)) {
    add <- lapply(changes, function(y) {y$prop_add})
    add <- add[!(is.na(add))]
    add <- Reduce(igraph::union, add)
    del <- lapply(changes, function(y) {y$prop_del})
    del <- del[!(is.na(del))]
    del <- Reduce(igraph::union, del)
  }
  if(is.null(del) & is.null(add)) {
    stop("Please supply modifications (additions and/or deletions")
  }
  g <- x$taxa

  # delete edges (this step is simple)
  if(length(del) > 0) {
    if(!is.numeric(del)) {del <- igraph::as_ids(del)}
    g <- igraph::delete.edges(g, edges = del)
  }
  # add edges (the new vertices must be made as a separate graph first to give directed edges)
  if(length(add) > 0) {
    add <- igraph::graph(rev(add$name))
    g <- g + add
  }
  # update vertex characteristics
  if(any(length(add) > 0 | length(del) > 0)) {
    igraph::V(g)$degree <- igraph::degree(g, mode = "out")
    igraph::V(g)$degree[which(igraph::V(g)$degree == 0)] <- 1
    igraph::V(g)$strength <- igraph::strength(graph = g, mode = "out", weights = igraph::E(g)$freq)
    igraph::V(g)$strength[which(igraph::V(g)$strength == 0)] <- 1
  }
  x$taxa <- g
  if(max(igraph::V(g)$degree) < 2) {
    x$status <- "Resolved"
  } else {
    x$status <- "Unresolved"
  }
  return(x)
}
