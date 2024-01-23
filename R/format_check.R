#' format_check
#'
#' Function to perform a series of basic formatting checks
#' geared towards taxonomic name data. The function very
#' simply checks for non letter characters in the taxonomic
#' names, that species-level names contain two words, and
#' genus-level and above names contain one word.
#' @param x A dataframe with hierarchically organised, taxonomic
#' information. If x only comprises the taxonomic information,
#' @param ranks does not need to be specified, but the columns
#' must be in order of decreasing taxonomic rank @param ranks
#' The column names of the taxonomic data fields in x. These
#' must be provided in order of decreasing taxonomic rank
#' @param species A logical indicating if x contains a species
#' column. As the data must be supplied in hierarchical order,
#' this column will naturally be the last column in x and
#' species-specific spell checks will be performed on this column.
#' @param species_sep A character vector of length one specifying
#' the genus name and specific epithet in the species column, if
#' present
#' @param verbose A logical determining if any flagged
#' errors should be reported to the console
#' @return A list of two lists. The first list flags the row
#' indexes of columns whose elements contains non-letter characters.
#' The second list flags the row indexes  of columns whose elements
#' do not contain the correct numbers of words
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # define ranks
#' b_ranks <- c("phylum", "class", "order", "family", "genus")
#' # run function
#' flag <- format_check(brachios, ranks = b_ranks)

format_check <- function (x, ranks, species = FALSE, species_sep = " ", verbose = TRUE) {

  # x = occ_data_raw_filtered
  # ranks = ranks
  # species = TRUE
  # species_sep = " "
  # verbose = TRUE

  if (!exists("x")) {
    stop("Please supply x as a dataframe of taxonomic assignments")
  }
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- as.character(1:ncol(x))
  }
  if (is.null(ranks)) {
    ranks <- colnames(x)
  }
  if (!all(ranks %in% colnames(x))) {
    stop("Not all elements of argument ranks are column names in x")
  }
  if (length(unique(x[, ranks[length(ranks)]])) < length(unique(x[,
                                                                  ranks[(length(ranks) - 1)]]))) {
    warning("Higher taxonomy is more diverse than lower taxonomy. Are the columns in x\n            or the column names specified in 'ranks' supplied in descending hierarchical order?")
  }
  if (!all(apply(x[, ranks], 2, class) == "character")) {
    stop("Not all columns in x are of class character")
  }
  if (!is.logical(species) & length(species) != 1) {
    stop("Species should be a logical of length one, indicating whether species-level designations are present in x")
  }
  if (species) {
    if (!is.character(species_sep) & length(species_sep) !=
        1) {
      stop("species_sep should be a character string identifying the genus and specific epithet separator in the species name column")
    }
  }
  x <- x[, ranks]
  chars <- list()
  chars2 <- vector()
  lens <- list()
  lens2 <- vector()
  if (species) {
    x[, length(ranks)] <- gsub(species_sep, " ", x[, length(ranks)])
    for (i in 1:length(ranks)) {
      chars[[i]] <- which(grepl("[^ [:alpha:]]", x[, ranks[i]]))
      chars2[i] <- ifelse(length(chars[[i]]) > 0, TRUE, FALSE)
    }
    for (i in 1:(length(ranks) - 1)) {
      lens[[i]] <- which(as.logical(unlist(lapply(strsplit(x[, ranks[i]], " "), length)) - 1))
      lens2[i] <- any(as.logical(unlist(lapply(strsplit(x[, ranks[i]], " "), length)) - 1))
    }
    lens[[length(ranks)]] <- which(as.logical(unlist(lapply(strsplit(x[,
                                                                       ranks[i]], " "), length)) - 2))
    lens2[length(ranks)] <- any(as.logical(unlist(lapply(strsplit(x[,
                                                                    ranks[i]], " "), length)) - 2))
    if (sum(chars2) != 0) {
      if (verbose) {
        message(paste0("Non-letter characters detected at the following ranks: ",
                       paste0(ranks[chars2], collapse = ", ")))
      }
    }
    if (sum(lens2[1:(length(ranks) - 1)]) != 0) {
      if (verbose) {
        message(paste0("The following ranks contain names consisting of more than one word: ",
                       paste0(ranks[lens2], collapse = ", "), ". Supraspecific taxon names are assumed to consist of single words"))
      }
    }
    if (lens2[length(ranks)] != 0) {
      if (verbose) {
        message(paste0("The species colum contain names consisting of more than two words. Species names are assumed to consist of two words"))
      }
    }
  } else {
    for (i in 1:length(ranks)) {
      chars[[i]] <- which(grepl("[^[:alpha:]]", x[, ranks[i]]))
      chars2[i] <- ifelse(length(chars[[i]]) > 0, TRUE, FALSE)
    }
    for (i in 1:length(ranks)) {
      lens[[i]] <- which(as.logical(unlist(lapply(strsplit(x[,
                                                             ranks[i]], " "), length)) - 1))
      lens2[i] <- any(as.logical(unlist(lapply(strsplit(x[,
                                                          ranks[i]], " "), length)) - 1))
    }
    if (sum(chars2) != 0) {
      if (verbose) {
        message(paste0("Non-letter characters detected at the following ranks: ",
                       paste0(ranks[chars2], collapse = ", ")))
      }
    }
    if (sum(lens2[1:(length(ranks) - 1)]) != 0) {
      if (verbose) {
        message(paste0("The following ranks contain names consisting of more than one word: ",
                       paste0(ranks[lens2], collapse = ", "), ". Supraspecific taxon names are assumed to consist of single words"))
      }
    }
  }
  names(chars) <- names(lens) <- ranks
  out <- list(chars, lens)
  names(out) <- c("non-letter", "word-count")
  return(out)
}
