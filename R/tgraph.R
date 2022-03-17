#' tgraph
#'
#' Function to create a tgraph representation of a
#' hierarchically organised dataframe. This is the focal
#' object of the t* functions - the complete set of
#' hierarchical relationships between a set of elements
#' @param x A dataframe containing a set of hierarchical
#' relationships. The leftmost column contains the
#' elements which will form the highest rank, followed
#' rightwards by successive ranks
#' @param ranks If not NULL, a vector of column names of x,
#' given in rank order. This is useful if x contains
#' columns which are not rank relevant or if columns are
#' not in hierarchical order. If not supplied, the
#' column order in x is used directly and is assumed to
#' be in rank order
#' @param verbose A logical indicating whether the
#' progress of tgraph construction should be reported
#' to the console
#' @return a tgraph object
#' @import igraph
#' @importFrom stats complete.cases

tgraph <- function(x, ranks = NULL, verbose = TRUE) {

  # check that the right object has been supplied
  if(!exists("x") | !(class(x)[1] %in% c("data.frame", "matrix", "array"))) {
    stop("Please provide a table of taxonomic classifications")
  }

  # get table of edges
  if(verbose) {
    cat(paste0("Filtering data"), "\n")
  }
  char <- apply(x, 2, is.character)
  if(any(isFALSE(char))) {
    stop("One or more of the columns contain non character data")
  }

  # check data for rank names (if ranks = NULL, the dataset colnames are used and are assumed to be in taxonomic order)
  if(is.null(ranks)) {
    ranks <- colnames(x)
    # otherwise use the ranks argument (assumed to be supplied in taxonomic order) is used to select the columns of interest
  } else {
    if(sum(ranks %in% colnames(x)) != length(ranks)) {
      stop("Not all specified ranks are present in the dataset")
    } else {
      x <- x[,ranks]
    }
  }

  # remove totally incomplete ranks
  ind <- apply(x, 1, function(x) all(is.na(x)))
  x <- x[!ind,]

  # flag non-discrete names if they occur
  test <- discrete_ranks(x = x, ranks = ranks)
  t1 <- unlist(test$crossed_adj)
  t2 <- unlist(test$crossed_all)
  if(length(t1) != 0 | length(t2) != 0) {
    stop("One or more elements of x are not unique to their rank")
  }

  # set the numeric definitions for each rank (1 is the largest rank in the dataset, smaller numbers give smaller ranks)
  rnum <- 1:length(ranks)
  names(rnum) <- ranks

  # record the number of columns (number of ranks)
  orig <- ncol(x)
  # create numeric ids for the rank of each column
  for(i in 1:ncol(x)) {
    x[, ncol(x) + 1] <- i
    names(x)[ncol(x)] <- paste0("L", i)
  }
  # for NA values, the value in the next highest column is used to ensure a complete structure across levels
  for(i in 2:orig) {
    x[,paste0("L", i)][which(is.na(x[,i]))] <- x[,paste0("L", i - 1)][which(is.na(x[,i]))]
    x[,i][which(is.na(x[,i]))] <- x[,(i - 1)][which(is.na(x[,i]))]
  }
  colnames(x) <- NULL

  # get table of edges
  if(verbose) {
    cat(paste0("Generating table of edges"), "\n")
  }
  gl <- list()
  for(i in 1:(orig - 1)) {
    gl[[i]] <- cbind(x[,c(c(i, i + 1), c(orig + i, orig + i + 1))])
  }
  g_list <- do.call(rbind, gl)
  # for taxa without links, a placeholder is added, for removal in the graph
  g_list[is.na(g_list[,1]), 3] <- NA
  g_list[is.na(g_list[,1]), 1] <- "ORPHAN"
  g_list[is.na(g_list[,2]), 4] <- NA
  g_list[is.na(g_list[,2]), 2] <- "ORPHAN"
  # remove loops
  lp <- which(g_list[,1] == g_list[,2])
  if(length(lp) > 0) {
    g_list <- g_list[-lp, ,drop = FALSE]
  }

  # get the frequencies of links before trimming
  if(verbose) {
    cat(paste0("Tabulating edge frequencies"), "\n")
  }
  g_val <- paste(g_list[,1], g_list[,2], sep = "|")
  g_val <- table(g_val)
  # get the unique links (g_list)
  g_list <- g_list[!duplicated(g_list[,1:2]),]
  g_val2 <- paste(g_list[,1], g_list[,2], sep = "|")
  # recover the frequencies of the unique links
  g_list$freq <- g_val[match(g_val2, names(g_val))]
  tax_c <- table(g_list[,4])
  g_list$q <- as.vector(tax_c[match(g_list[,4], names(tax_c))])
  # get the name datatable
  g_name <- as.data.frame(cbind(c(g_list[,1], g_list[,2]),
                                c(g_list[,3], g_list[,4])))
  g_name <- unique(g_name)
  # make the graph (convenient way of calculating taxon statistics)
  if(verbose) {
    cat(paste0("Building tgraph"), "\n")
  }
  t_graph <- igraph::graph_from_edgelist(as.matrix(g_list[,2:1]))
  # edge weights (inverse used for pathways so that more common edges have lower weights i.e. less cost)
  # normal value (freq) used to permit calculation of vertex strength and its derived measure, the highest edge strength proportion
  igraph::E(t_graph)$weight <- 1 / g_list[,5]
  igraph::E(t_graph)$freq <- g_list[,5]
  # vertex character and numeric ranks
  igraph::V(t_graph)$rank <- as.numeric(g_name[(match(igraph::V(t_graph)$name, g_name[,1])),2])
  igraph::V(t_graph)$inf_rank <- NA
  # remove the placeholder ORPHAN vertex
  if("ORPHAN" %in% g_list[,1] | "ORPHAN" %in% g_list[,2]) {
    t_graph <- igraph::delete.vertices(t_graph, "ORPHAN")
  }
  g_list <- g_list[-(which(g_list[,1] == "ORPHAN" | g_list[,2] == "ORPHAN")),]
  g_name <- g_name[stats::complete.cases(g_name),]
  # number of vertices that the focal vertex links out to
  igraph::V(t_graph)$degree <- igraph::degree(t_graph, mode = "out")
  # vertex strength (number of times that vertex links out in the database, 0 also assigned to 1)
  igraph::V(t_graph)$strength <- igraph::strength(graph = t_graph, mode = "out", weights = igraph::E(t_graph)$freq)

  # output tgraph object
  out <- list()
  out[[1]] <- t_graph
  out[[2]] <- rnum
  out[[3]] <- "graph"
  if(max(unique(igraph::V(t_graph)$degree)) > 1) {
    out[[4]] <- "Unresolved"
  } else {
    out[[4]] <- "Resolved"
  }
  names(out) <- c("taxa", "ranks", "type", "status")
  class(out) <- "tgraph"
  return(out)
}
