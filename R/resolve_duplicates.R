#' resolve_duplicates
#'
#' Function for identifying and resolving alternative higher
#' assignments in a hierarchically structured dataframe.
#' Columns are checked from the lowest to the highest rank for
#' elements with multiple higher assignments. These assignments
#' are then assessed topologically to determine if they
#' represent inadvertent use of the same name at a given rank
#' for genuinely different entities, or whether the higher
#' classifications are conflicting. In the case of the former,
#' unique character suffixes are applied to each differently
#' classified case (up to 26 currently supported), effectively
#' splitting up the alternatively classified element. In the
#' case of the latter, the alternative classifications are
#' assessed and are either combined, or the more frequently
#' used or the more complete classification scheme is taken
#' (the more frequent pathway can also be the most complete).
#' @param x A dataframe containing hierarchically structured
#' information, for example a table of genus names and their
#' higher taxonomic classifications
#' @param ranks If not NULL, a vector of column names of x,
#' given in rank order. This is useful if x contains columns
#' which are not rank relevant or if columns are not in
#' hierarchical order. If not supplied, the column order in x
#' is used directly and is assumed to be in rank order
#' @param jump The maximum number of levels between the point of
#' divergence and the point of reunion (if present) for a given
#' path, below which the divergence will be taken as conflicting
#' @param plot A logical speciying if the divergent paths should
#' be plotted
#' @param verbose A logical of length one which determines if
#' the function should report the detection and resolution of
#' elements with multiple higher classifications (if any)
#' @return The dataframe x, with any alternative higher
#' classifications resolved, giving the classification a strict
#' tree structure
#' @import igraph
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # define ranks
#' b_ranks <- c("phylum", "class", "order", "family", "genus")
#' # run function
#' res <- resolve_duplicates(brachios, ranks = b_ranks)

resolve_duplicates <- function(x, ranks = NULL, jump = 4, plot = FALSE, verbose = TRUE) {

  # check arguments
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
  if(!all(unlist(lapply(x[,levels], class)) == "character")) {
    stop("Not all rank columns in x are of class character")
  }

  for(i in length(ranks):2) {

    if(i != length(ranks)) {cat("\r")}
    cat(paste0(" + resolving duplicates at rank ", ranks[i], "      "))
    if(i == 2) {cat("\n")}

    link <- unique(x[,ranks[c(1:i)]])
    tcount <- table(link[,ncol(link)])
    whichd <- names(tcount)[tcount > 1]
    if(length(whichd) > 0) {

      # notify R if specified
      if(verbose) {
        #cat(sprintf("%d duplicates detected at rank %s. Resolving...", length(whichd), ranks[i]), "\n")
      }
      for(j in 1:length(whichd)) {

        # get all instances of the duplicate, along with their higher classifications
        subst <- x[which(x[,ranks[i]] == whichd[j]), ranks[1:which(ranks == ranks[i])]]
        # small tgraph
        st <- tgraph(subst, verbose = FALSE)
        # coerce to tvertseq (mode = parent)
        sv <- list()
        sv[[1]] <- igraph::induced_subgraph(st$taxa, which(igraph::V(st$taxa)$rank == max(st$ranks)))
        sv[[1]] <- igraph::delete.edges(sv[[1]], edges = igraph::E(sv[[1]]))
        sv[[2]] <- st$ranks
        sv[[3]] <- st$taxa
        sv[[4]] <- mode
        names(sv) <- c("taxa", "ranks", "seq", "mode")
        class(sv) <- "tvertseq"

        # if the 'duplicates' are cases where the higher ranks are NA
        if(igraph::V(sv$taxa)$degree == 1) {

          # derive the complete higher taxonomy
          higher_v <- igraph::ego(st$taxa, nodes = whichd[j], mode = "out",
                          order = max(igraph::V(st$taxa)$rank))[[1]]
          higher <- names(higher_v)
          higher_r <- igraph::V(st$taxa)[higher_v]$rank
          highf <- rep(NA, max(igraph::V(st$taxa)$rank))
          highf[higher_r] <- higher
          x[which(x[,ranks[i]] == whichd[j]), ranks[1:max(higher_r)]] <-
            rep(highf, each = length(which(x[,ranks[i]] == whichd[j])))

          # if the focal taxa has true duplicate classifications
        } else {
          # resolve duplicate paths
          to_do <- assess_duplicates(st, node = whichd[j], jump = jump)

          # either update taxonomy
          if(length(grep("proposed", to_do[[1]]$dec)) == 1) {
            to_do <- update_graph(st, changes = to_do)
            # get the immediate higher classification
            higher_v <- igraph::ego(st$taxa, nodes = whichd[j], mode = "out",
                            order = max(igraph::V(st$taxa)$rank))[[1]]
            higher <- names(higher_v)
            higher_r <- igraph::V(st$taxa)[higher_v]$rank
            highf <- rep(NA, max(igraph::V(st$taxa)$rank))
            highf[higher_r] <- higher
            x[which(x[,ranks[i]] == whichd[j]), ranks[1:max(higher_r)]] <-
              rep(highf, each = length(which(x[,ranks[i]] == whichd[j])))

            # or add suffixes for each unique classification
          } else {
            sbr <- unique(subst)
            # remove all NA combinations
            nat <- apply(sbr, 1, function(x) {sum(is.na(x))})
            rm_nat <- which(nat == (ncol(sbr) - 1))
            if(length(rm_nat) == 1) {
              sbr <- sbr[-rm_nat,]
            }
            for(k in 1:nrow(sbr)) {
              # index last non NA higher classification. match name at that rank
              pos <- max(which(!is.na(sbr[k,1:(ncol(sbr) - 1)])))
              x[which(x[,colnames(sbr)[pos]] == sbr[k,pos] &
                        x[,ranks[i]] == sbr[1,ncol(sbr)]), ranks[i]] <- paste0(whichd[j], LETTERS[k])
            }
          }
        }
      }
    }
  }
  return(x)
}
