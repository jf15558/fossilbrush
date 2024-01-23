#' roxygen documentation
#'
#' discrete_ranks
#'
#' Function for checking whether names in one column of a
#' hierarchically organised dataframe re-occur at other levels.
#' Two checks are performed. The first checks for names in
#' adjacent column, assuming that accidental reuse of names at
#' other levels are most likely to occur at an adjacent rank.
#' The second compares across all columns.
#'
#' @param x A dataframe containing hierarchically structured
#' information, for example a table of genus names and their
#' higher taxonomic classifications
#' @param ranks If not NULL, a vector of column names of x,
#' given in rank order. This is useful if x contains columns
#' which are not rank relevant or if columns are not in
#' hierarchical order. If not supplied, the column order in x
#' is used directly and is assumed to be in rank order
#' @return A list of two lists. The first list contains names
#' which reoccur at adjacent ranks. The second list contains
#' names that reoccur at any rank
#' @importFrom stats na.omit
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # define ranks
#' b_ranks <- c("phylum", "class", "order", "family", "genus")
#' # run function
#' flag <- discrete_ranks(brachios, ranks = b_ranks)

discrete_ranks <- function(x, ranks = NULL) {

  if(!exists("x")) {
    stop("Please supply a dataframe and optionally a vector of ranks to check")
  }
  if(!is.data.frame(x)) {
    stop("x must be a dataframe")
  }
  if(!is.null(ranks)) {
    if(!all(ranks %in% colnames(x))) {
      stop("Not all elements in ranks are column names of x")
    }
    levels <- ranks
  } else {
    levels <- colnames(x)
  }
  x <- x[,levels]
  if(!all(unlist(lapply(x, class)) == "character")) {
    stop("Not all rank columns in x are of class character")
  }

  crossed_adj <- list()
  crossed_all <- list()
  for(i in 1:(length(levels) - 1)) {
    # get the cross level taxa
    crossed_adj[[i]] <- unique(na.omit(x[,levels[i]]))[unique(na.omit(x[,levels[i]])) %in% unique(na.omit(x[,levels[i + 1]]))]
    crossed_all[[i]] <- unique(na.omit(x[,levels[i]]))[unique(na.omit(x[,levels[i]])) %in% unique(na.omit(as.vector(unlist(x[,levels[-i]]))))]
  }
  names(crossed_adj) <- paste0(levels[-length(levels)], "--", levels[-1])
  names(crossed_all) <- levels[-length(levels)]
  out <- list(crossed_adj, crossed_all)
  names(out) <- c("crossed_adj", "crossed_all")
  return(out)
}
