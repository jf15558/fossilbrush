#' plot_taxa
#'
#' Function to plot the parent or child relationships of an
#' element in a hierarchically organised dataframe. Multiple
#' taxa can be plotted simultaneously
#' @param x a dataframe containing hierarchically organised
#' data in columns
#' @param taxon A character vector of element names whose
#' relationships will be plotted (these must be of the same
#' rank)
#' @param trank A character vector of length one corresponding
#' to the column name in x in which taxa is located
#' @param ranks A character vector corresponding to the column
#' names in x, given in hierarchical order
#' @param mode The direction of the relationships to be
#' plotted
#' @param step A positive integer specifinyg the
#' neighbourhood of the relationships to plot. Specifying
#' a number greater than the number of ranks will not cause a
#' failure, and will instead plot all relationships in the
#' direction specified in mode
#' @return A plot of the relationships of the specified
#' elements
#' @import igraph
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # define ranks in dataset
#' b_ranks <- c("phylum", "class", "order", "family", "genus")
#' # plot taxon
#' plot_taxa(brachios, "Atrypa", trank = "genus", ranks = b_ranks, mode = "parent")

plot_taxa <- function(x, taxon, trank, ranks, mode = c("parent", "child", "all"), step = NULL) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # check that data has minimally been supplied
  if(!exists("x")) {
    stop("Please supply x as a dataframe of taxonomic assignments")
  }
  if(!is.data.frame(x)) {
    stop("x must be a dataframe. For plotting t* objects directly, use plot_t")
  }
  # check taxon data
  if(!exists("taxon")) {
    stop("Please provide at least one taxon to plot relationships for")
  }
  if(!is.character(taxon)) {
    stop("taxon should be a vector of mode character")
  }
  # check rank data
  if(!exists("trank")) {
    stop("Please provide the rank(s) of the taxa to be plotted")
  }
  if(!is.character(trank) | length(trank) != length(taxon)) {
    stop("trank should be a vector of mode character the same length as taxon")
  }
  if(!exists("ranks")) {
    stop("Please provide ranks as the taxonomic column names in x")
  }
  if(!is.character(ranks)) {
    stop("ranks should be a vector of mode character")
  }
  if(!all(c(trank, ranks) %in% colnames(x))) {
    stop("all elements of trank and ranks must be column names in x")
  }
  if(!all(apply(x[,ranks], 2, class) == "character")) {
    stop("Not all columns in x are of class character")
  }
  if(length(mode) != 1) {
    stop("Please specify the plotting mode as one of: all, parent, child")
  }
  if(!mode %in% c("all", "parent", "child")) {
    stop("mode must be one: all, parent or child")
  }
  if(mode == "parent") {
    mode <- "out"
  }
  if(mode == "child") {
    mode <- "in"
  }
  if(!is.null(step)) {
    if(!is.numeric(step) | length(step) != 1) {
      stop("step should be a numeric of length 1 specifying the neighbourhood of the relationships to plot")
    }
  } else {
    step <- length(ranks)
  }

  # subset to columns
  x <- x[,ranks]
  # check that the data is character
  if(!all(apply(x[,ranks], 2, class) == "character")) {
    stop("Not all columns in x are of class character")
  }

  # tgraph from the subset of the dataframe containing the focal taxa
  tlist <- list()
  for(i in 1:length(taxon)) {
    if(!taxon[i] %in% x[,trank[i]]) {
      warning(paste0("element ", i, " in taxon is not present the column of x denoted by its corresponding element in trank and will not be plotted"))
    }
    tlist[[i]] <- unique(x[which(x[,trank[i]] == taxon[i]), ranks])
  }
  tg <- tgraph(do.call(rbind, tlist), verbose = FALSE)
  tnode <- which(igraph::V(tg$taxa)$name %in% taxon)

  try <- igraph::ego(tg$taxa, order = step, nodes = tnode, mode = mode)
  try <- lapply(try, function(y) {igraph::induced_subgraph(tg$taxa, y)})
  try <- igraph::V(Reduce(igraph::union, try))$name
  try <- igraph::induced_subgraph(tg$taxa, try)

  labs <- rep("", length(unique(igraph::V(try)$rank)))
  labs[as.numeric(as.factor(igraph::V(try)$rank))] <- names(tg$ranks)[igraph::V(try)$rank]
  rseq <- seq(from = -1, to = 1, length.out = length(labs))

  igraph::V(try)$color = rep(igraph::categorical_pal(8)[1], length(igraph::V(try)))
  igraph::V(try)$rank <- as.numeric(as.factor(igraph::V(try)$rank))
  igraph::V(try)$color[igraph::V(try)$name %in% taxon] <- igraph::categorical_pal(8)[2]
  lay <- igraph::layout_with_sugiyama(try, layers = as.numeric(igraph::V(try)$rank), hgap = 1)

  par(mar = c(0, 5, 0, 5))
  par(las = 1)
  plot(try, layout = lay$layout, edge.arrow.mode = 0, margin = 0.2)
  axis(4, labels = rev(labs), at = rseq, col = NA)
}
