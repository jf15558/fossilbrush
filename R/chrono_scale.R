#' chrono_scale
#'
#' Convenience function to apply user-specified chronostratigraphy
#' to fossil datasets. The function relies on a lookup table
#' generated based on the named intervals in the PBDB in early
#' 2021. First and last interval names in the supplied dataset
#' are matched against this lookup table, by default using
#' 'get("GTS_2020)", to get GTS_2020 numeric ages. If the
#' dataset contains intervals which are not present in the
#' lookup table, they will not be matched and the user will
#' be warned. To get around this possibility, the user can
#' also supply the original numeric ages which will be used
#' as default ages if an interval cannot be matched, to ensure
#' that the returned vectors of numeric ages do not contain NAs.
#' @param x A data.frame containing, minimally two columns
#' corresponding respectively to the first and last intervals
#' of the data. Values should only be present in the second column
#' where the minimum age interval for a row is different to the
#' maximum age interval. Otherwise the values should be NA and the
#' ages returned will be based on the interval specified in the
#' first column, in line with PBDB formatting.
#' @param tscale A character string specifying one of the inbuilt
#' chronostratigraphic timescales (currently GTS 2020 only) or a
#' data.frame supplied by the user. If the latter, this must
#' contain columns named 'Interval', 'FAD', 'LAD', specifying the
#' interval names to be matched and their lower and upper age in Ma
#' @param srt A character of length 1 specifing the column name
#' of the first interval field in x
#' @param end A character of length 1 specifing the column name
#' of the last interval field in x
#' @param max_ma If not NULL, a character of length 1 specifing
#' the column name of the original numeric maximum age field in
#' x, to be used as fall back values if interval names cannot
#' all be matched
#' @param min_ma If not NULL, a character of length 1 specifing
#' the column name of the original numeric minimum age field in
#' x, to be used as fall back values if interval names cannot
#' all be matched
#' @param verbose A logical indicating if warning messages should be
#' displayed or otherwise
#' @return The dataframe, x, with two additional columns containing
#' the revised first and last numeric ages of the data, with column
#' names GTS_FAD and GTS_LAD respectively
#' @export
#' @examples
#' # example dataset
#' data("brachios")
#' # add GTS_2020 dates
#' brachios <- chrono_scale(brachios, srt = "early_interval", end = "late_interval",
#'                           max_ma = "max_ma", min_ma = "min_ma")

chrono_scale <- function(x, tscale = "GTS_2020", srt = "early_interval", end = "late_interval", max_ma = NULL, min_ma = NULL, verbose = TRUE) {

  if(!exists("x")) {
    stop("Please supply x as a dataframe containing, minimally, the first and last interval ages of PBDB data")
  }
  if(!is.data.frame(x)) {
    stop("Please supply x as a dataframe containing, minimally, the first and last interval ages of PBDB data")
  }
  if(is.character(tscale)) {
    if(length(tscale) != 1) {
      stop("Only a single chronostratigraphic dataset name should be supplied")
    }
    if(!tscale %in% c("GTS_2020")) {
      stop("tscale is not a valid chronostratigraphic dataset name")
    }
    if(tscale == "GTS_2020") {tscale <- fossilbrush::GTS_2020}

  } else if(is.data.frame(tscale)) {
    if(ncol(tscale) < 3) {
      stop("If supplying a dataframe, this must contain Interval, FAD and LAD columns")
    }
    if(!all(c("Interval", "FAD", "LAD") %in% colnames(tscale))) {
      stop("If supplying a dataframe, this must contain Interval, FAD and LAD columns")
    }
    tscale <- tscale[, c("Interval", "FAD", "LAD")]
    if(any(is.na(tscale))) {
      stop("NA values are not permitted in tscale, if supplying a data.frame")
    }
    if(!is.character(tscale[,"Interval"])) {
      stop("The interval column in tscale should contain interval names")
    }
    if(!all(c(is.numeric(tscale[,"FAD"]), is.numeric(tscale[,"LAD"])))) {
      stop("The FAD and LAD columns in tscale should be numeric")
    }
    if(any(tscale[,c("FAD", "LAD")] < 0)) {
      stop("One or more entries in the FAD or LAD columns is negative")
    }
    if(any(tscale[,"FAD"] - tscale[,"LAD"] <= 0)) {
      stop("One or more LADs are the same age or older than their corresponding FADs")
    }

  } else {
    stop("Supply tscale as a valid chronostratigraphic dataset name or a data.frame")
  }
  if(length(srt) != 1 | !is.character(srt)) {
    stop("srt should be a character vector of length 1 and refer to the first interval ages in x")
  }
  if(!srt %in% colnames(x)) {
    stop("The supplied value of srt is not a column name in x")
  }
  if(!is.character(x[,srt])) {
    stop("The value of srt does not refer to a character class column in x")
  }
  if(length(end) != 1 | !is.character(end)) {
    stop("end should be a character vector of length 1 and refer to the last interval ages in x")
  }
  if(!end %in% colnames(x)) {
    stop("The supplied value of end is not a column name in x")
  }
  if(!all(is.na(x[,end]))) {
    if(!is.character(x[,end])) {
      stop("The value of end does not refer to a character class column in x")
    }
  }
  if(is.null(max_ma) != is.null(min_ma)) {
    stop("If supplying the original ages, both max_ma and min_ma must be specified")
  }
  if(!is.null(max_ma)) {

    if(length(max_ma) != 1 | !is.character(max_ma)) {
      stop("max_ma should be a character vector of length 1 and refer to the first interval ages in x")
    }
    if(!max_ma %in% colnames(x)) {
      stop("The supplied value of max_ma is not a column name in x")
    }
    if(!is.numeric(x[,max_ma])) {
      stop("The value of max_ma does not refer to a numeric class column in x")
    }
    if(length(min_ma) != 1 | !is.character(min_ma)) {
      stop("min_ma should be a character vector of length 1 and refer to the last interval ages in x")
    }
    if(!min_ma %in% colnames(x)) {
      stop("The supplied value of min_ma is not a column name in x")
    }
    if(!is.numeric(x[,min_ma])) {
      stop("The value of min_ma does not refer to a numeric class column in x")
    }
  }

  # set up vectors
  xfad <- x[,srt]
  xlad <- x[,end]
  xlad[is.na(xlad)] <- xfad[is.na(xlad)]
  xerl <- x[,max_ma]
  xlte <- x[,min_ma]
  cinterval <- "Interval"
  cfad <- "FAD"
  clad <- "LAD"
  # check for unmatchable names
  test <- unique(c(xfad, xlad))
  test <- test[!test %in% tscale[,cinterval]]
  if(length(test) > 0) {
    if(verbose) {warning(paste0("The following intervals are not present in the GTS_2020 lookup: ", paste0(test, collapse = ", "), ". If max_ma and min_ma have been supplied, these values will be supplied for the unmatched intervals as defaults. Otherwise, NA will be returned for these intervals"))}
  }

  new_fad <- tscale[match(xfad, tscale[,cinterval]), cfad]
  new_lad <- tscale[match(xlad, tscale[,cinterval]), clad]
  if(!is.null(max_ma)) {
    new_fad[is.na(new_fad)] <- xerl[is.na(new_fad)]
    new_lad[is.na(new_lad)] <- xlte[is.na(new_lad)]
  }

  # ensure that the returned ages are the right way around
  news <- cbind(new_fad, new_lad)
  newFad <- apply(news, 1, function(x) {if(all(is.na(x))) {return(NA)} else {return(which.max(x))}})
  newLad <- apply(news, 1, function(x) {if(all(is.na(x))) {return(NA)} else {return(which.min(x))}})

  x$newFAD <- newFad
  x$newLAD <- newLad
  return(x)
}
