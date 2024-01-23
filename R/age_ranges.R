#' roxygen documentation
#'
#' age_ranges
#'
#' Function to derive a range table of taxon names from a
#' stratigraphic occurrence dataset. The default behaviour
#' is to return a total range table - the oldest FAD and
#' youngest LAD for each taxon (max), but the function can
#' also return the minimum range - youngest FAD and oldest
#' LAD (min), or the uncertainty bounds on each FAD and
#' LAD - the two oldest FADs and two youngest LADs (bounds).
#' The names for which ranges are derived are specified by
#' the taxononmy argument, but multiple elements can be
#' given here, allowing taxonomic range for higher clades
#' to also be returned.
#' @param data A three column dataframe comprising one or
#' more character columns of taxonomic names, a numeric
#' column of FADs and a numeric column of LADs
#' @param taxonomy A character vector corresponding to one or
#' more of the taxonomic name columns in data
#' @param srt A character vector of length one specifying the
#' FAD column in data
#' @param end A character vector of length one specifying the
#' LAD column in data
#' @param mode A character vector of length one specifying
#' the type of range table to return: one of max, min or
#' bounds. If not specified by the user, the function
#' behaviour will default to max
#' @return A dataframe containing at least four columns:
#' taxon name, FAD, LAD and the taxonomic rank. If taxonomy
#' is of length one, taxonomic rank will be a vector of
#' identical names. If mode = "bounds", there will be two
#' pairs of age columns, denoting the upper and lower bounds
#' on the FAD and LAD for each taxon name
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # derive age ranges
#' rng <- age_ranges(brachios)

age_ranges <- function(data, taxonomy = "genus", srt = "max_ma", end = "min_ma", mode = "max") {

  # check arguments
  if(!is.data.frame(data)) {
    stop("Data should be a dataframe")
  }
  if(!all(c(taxonomy, srt, end) %in% colnames(data))) {
    stop("One or more of taxonomy, srt or end are not colnames in data")
  }
  if(!all(mode %in% c("max", "min", "bounds"))) {
    stop("Mode must be one of the following: max, min, bounds")
  }
  if(length(mode) > 1) {
    stop("Mode must be one of the following: max, min, bounds")
  }
  if(!is.numeric(data[,srt]) | !is.numeric(data[,end])) {
    stop("Columns srt and end must be numeric")
  }
  if(any(is.na(data[,srt])) | any(is.na(data[,end]))) {
    stop("One or more of the ages is NA")
  }
  if(any(data[,srt] < data[,end])) {
    stop("One or more maximum ages are smaller than their corresponding minimum ages")
  }

  # for a single rank
  if(length(taxonomy) == 1) {

    # add a small constant to the FAD to prevent errors in zero-range ages (i.e. FAD == LAD)
    data[,srt] <- data[,srt] + 0.1
    # reformat for indexing
    data <- cbind.data.frame(c(data[,taxonomy], data[,taxonomy]), c(data[,srt], data[,end]), c(rep("s", nrow(data)), rep("e", nrow(data))))
    data <- unique(data[complete.cases(data),])
    data <- data[order(data[,1], data[,2], decreasing = c(FALSE, TRUE), method = "radix"),]
    # unique positions (FAD uncertainty will be the two oldest unique FADs, same principle with LADS), achieved by ordering and indexing
    fad <- which(!duplicated(data[,1]))
    names(fad) <- data[fad,1]
    lad <- table(data[,1])
    lad <- lad[order(match(names(lad), names(fad)))]
    lad <- cumsum(lad)
    # remove constant
    data[which(data[,3] == "s"),2] <- data[which(data[,3] == "s"),2] - 0.1
    ages <- cbind.data.frame(data[fad,1], data[fad,2], data[(fad + 1),2], data[(lad - 1),2], data[lad,2], taxonomy)
    colnames(ages) <- c("taxon", "FAD_early", "FAD_late", "LAD_early", "LAD_late", "level")

    if(mode == "min") {
      ages <- ages[,c("taxon", "FAD_late", "LAD_early", "level")]
      colnames(ages) <- c("taxon", "FAD", "LAD", "level")
    }
    if(mode == "max") {
      ages <- ages[,c("taxon", "FAD_early", "LAD_late", "level")]
      colnames(ages) <- c("taxon", "FAD", "LAD", "level")
    }

    # for multiple ranks
  } else {
    ages <- list()
    for(i in 1:length(taxonomy)) {
      ages[[i]] <- age_ranges(data = data, taxonomy = taxonomy[i], srt = srt, end = end, mode = mode)
    }
    ages <- do.call(rbind.data.frame, ages)
  }
  return(ages)
}
