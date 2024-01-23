#' add_kingdoms
#'
#' Convenience function to add in a kingdom column to a PBDB
#' dataset. This relies on the dataset having a column of
#' phylum-level assignments for occurrences. The kingdom
#' column is a useful addition for filtering very large
#' taxonomically diverse datasets, and adds in an additional
#' level of data which can inform taxonomic cleaning routines
#' like those called by @seealso check_taxonomy
#' @param x A dataframe containing, minimally, phylum-level
#' assignments of the data
#' @param phylum A character of length 1 specifing the column
#' in x with the phylum level assignments
#' @param insert.left A convenience argument which will
#' make sure that the kingdom column will be inserted in
#' dataframe left immediately to the left of the phylum column
#' @return The dataframe x, with the kingdom column inserted
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # add kingdoms to dataset
#' brachios <- add_kingdoms(brachios)

add_kingdoms <- function(x, phylum = "phylum", insert.left = TRUE) {

  if(!exists("x")) {
    stop("Please supply x as a dataframe containing, minimally, the phylum-level assignments of PBDB data")
  }
  if(!is.data.frame(x)) {
    stop("Please supply x as a dataframe containing, minimally, the phylum-level assignments of PBDB data")
  }
  if(!is.character(phylum) | length(phylum) != 1) {
    stop("Phylum should a character vector of length 1 and refer to the phylum column in x")
  }
  if(!phylum %in% colnames(x)) {
    stop("The supplied value of phylum is not a column name in x")
  }
  if(!is.character(x[,phylum])) {
    stop("The column in x referred to by phylum is not of class character")
  }
  if("kingdom" %in% colnames(x)) {
    stop("x already contains a column called kingdoms")
  }
  if(!is.logical(insert.left)) {
    stop("insert.left should be a logical of length 1")
  }

  pbdb_kingdoms <- get("pbdb_kingdoms")
  animals <- pbdb_kingdoms$animals
  plants <- pbdb_kingdoms$plants
  protists <- pbdb_kingdoms$protists
  x$kingdom <- NA
  x$kingdom[x[,phylum] %in% animals] <- "Animalia"
  x$kingdom[x[,phylum] %in% plants] <- "Plantae"
  x$kingdom[x[,phylum] %in% protists] <- "Protista"
  x$kingdom[x[,phylum] == "Pezizomycotina"] <- "Fungi"
  x$kingdom[x[,phylum] == "Cyanobacteria"] <- "Bacteria"
  if(insert.left) {
    if(which(colnames(x) == "phylum") == 1) {
      x <- cbind.data.frame(kingdom = x$kingdom, x[,-which(colnames(x) == "kingdom"), drop = FALSE])

    } else {
      x <- cbind.data.frame(x[,1:(which(colnames(x) == phylum) - 1)], x$kingdom, x[,(which(colnames(x) == phylum)):(ncol(x) - 1)])
    }
  }
  colnames(x)[which(colnames(x) == "x$kingdom")] <- "kingdom"
  return(x)
}
