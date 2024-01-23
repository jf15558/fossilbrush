#' find_duplicates
#'
#' Function to detect and report elements with multiple
#' higher assigments in a hierarchically structured
#' dataframe
#' @param x A hierarchically organised dataframe
#' @param ranks The ranks in the dataframe in which
#' to check for elements with multiple higher
#' classifications. The top rank is ignored by default
#' @return A dataframe of elements with multiple
#' higher classifications and their ranks
#' @importFrom stats complete.cases
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' b_ranks <- c("phylum", "class", "order", "family", "genus")
#' # run function
#' flag <- find_duplicates(brachios, ranks = b_ranks)

find_duplicates <- function(x, ranks = NULL) {

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
  out <- data.frame(NA, NA)
  colnames(out) <- c("taxon", "rank")
  for(i in length(ranks):2) {
    whichd2 <- data.frame(NA, NA)
    colnames(whichd2) <- c("taxon", "rank")
    link <- unique(x[,ranks[c(1:i)]])
    tcount <- table(link[,ncol(link)])
    whichd <- names(tcount)[tcount > 1]
    if(length(whichd) > 0) {
      whichd2 <- cbind(whichd, rep(ranks[i], times = length(whichd)))
      colnames(whichd2) <- c("taxon", "rank")
    }
    out <- rbind.data.frame(out, whichd2)
  }
  out <- out[complete.cases(out),]
  return(out)
}
