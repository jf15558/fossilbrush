#' flag_ranges
#'
#' Function to compare stratigraphic ranges in x to a
#' set of reference ranges from y. A list of two elements
#' is returned. The first is a dataframe summarising
#' the overall error status, specific error counts
#' FAD and LAD differences, and the 95% density
#' distributions of the FAD and LAD errors for each
#' unique taxon in the column of x denoted by the first
#' element of xcols. If a taxon in x is not present in
#' y, it is assigned the status 000 and its other
#' entries in the returned dataframe will be NA.
#' The second element of the returned list is the error
#' code for every individual element of the column of
#' x denoted by the first element of xcols - this will
#' have the same number of rows as x. If x is a range
#' table rather than an occurrence dataset, then the two
#' list elements will have the same number of rows.
#' Ranges for comparison may be supplied directly in y,
#' or y may be another occurrence dataset, in which case
#' @seealso age_ranges is called internally to generate
#' the range table for comparison.
#' @param x Stratigraphic range data for taxa as a whole
#' or for individual fossil occurrences
#' @param y The same as in x. This is the dataset to
#' which ranges will be compared
#' @param xcols A character vector of length three
#' specifying, in the following order, the taxonomic name,
#' stratigraphic base (FAD) and stratigraphic top (LAD)
#' columns in x.
#' @param ycols An optional character vector of length
#' three for the same column types as in xcols, but for
#' dataset y. This is useful if the column names differ
#' between the datasets
#' @param flag.diff A vector of thresholds, given in millions
#' of years which will be used to flag discrepancies
#' between occurrence FADs and LADs with respect to the
#' reference range. This is a convenience parameter so that
#' occurrences with large discrepancies can be quickly
#' identified. Multiple thresholds can be supplied
#' @param verbose A logical of length one determining if
#' the flagging progress should be reported to the console
#' @return A list of two data.frames, the first recording
#' overall error statistics, the second recording error
#' types for each element of x. In the second data.frame,
#' FAD or LAD differences in excess of the supplied threshold(s)
#' are marked with 1, otherwise 0
#' @importFrom stats density
#' @import data.table
#' @export
#' @examples
#' # load the example datasets
#' data(brachios)
#' data(sepkoski)
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # update brachios to GTS2020 to match Sepkoski
#' brachios <- chrono_scale(brachios, srt = "early_interval", end = "late_interval",
#'                           max_ma = "max_ma", min_ma = "min_ma", verbose = FALSE)
#' brachios$max_ma <- brachios$newFAD
#' brachios$min_ma <- brachios$newLAD
#' # drop occurrences with older LADs than FADs
#' brachios <- brachios[brachios$max_ma > brachios$min_ma,]
#' # trim the Sepkoski Compendium to the relevant entries
#' sepkoski <- sepkoski[which(sepkoski$PHYLUM == "Brachiopoda"),]
#' # run flag ranges
#' flg <- flag_ranges(x = brachios, y = sepkoski, ycols = c("GENUS", "RANGE_BASE", "RANGE_TOP"))

flag_ranges <- function(x = NULL, y = NULL, xcols = c("genus", "max_ma", "min_ma"), ycols = NULL, flag.diff = 5, verbose = TRUE) {

  # check all information is supplied
  if(is.null(x) | is.null(y)) {
    stop("Please supply both x and y arguments")
  }
  if(!is.data.frame(x) | !is.data.frame(y)) {
    stop("x and y must both be dataframes")
  }
  # if the dataframes have less than three columns supplied, break
  if(ncol(x) < 3 | ncol(y) < 3) {
    stop("Supplied dataframes do not contain enough data (at least name, FAD and LAD needed in both")
  }
  if(length(xcols) != 3) {
    stop("xcols must contain three elements")
  }
  if(is.null(ycols)) {
    ycols <- xcols
  } else {
    if(length(ycols) != 3) {
      stop("ycols must contain three elements")
    }
  }
  if(!all(xcols %in% colnames(x))) {
    stop("One or more elements of xcols are not column names in x")
  }
  if(!all(ycols %in% colnames(y))) {
    stop("One or more elements of ycols are not column names in y")
  }
  if(!is.numeric(x[,xcols[2]]) | !is.numeric(x[,xcols[3]])) {
    stop("Elements 2 and 3 of xcols must refer to numeric columns in x")
  }
  if(!is.numeric(y[,ycols[2]]) | !is.numeric(y[,ycols[3]])) {
    stop("Elements 2 and 3 of ycols must refer to numeric columns in y")
  }
  if(!all(length(verbose) == 1, is.logical(verbose))) {
    stop("verbose must be one of TRUE or FALSE")
  }
  # global variable workaround for data.table syntax
  . <- NULL

  # subselect if many columns are present
  if(ncol(x) > 3) {
    x <- x[,xcols]
  }
  if(ncol(y) > 3) {
    y <- y[,ycols]
  }
  colnames(x) <- colnames(y) <- c("taxon", "max", "min")

  # count unique taxon names
  n1 <- unique(x[,"taxon"])
  n2 <- unique(y[,"taxon"])

  # derive checking range table if needed
  if(length(n1) != nrow(x)) {
    xr1 <- age_ranges(data = x, taxonomy = "taxon", srt = "max", end = "min")
  } else {
    xr1 <- x
  }
  xr1 <- xr1[order(xr1[,1]),]

  if(length(n2) != nrow(y)) {
    yr <- age_ranges(data = y, taxonomy = "taxon", srt = "max", end = "min")
  } else {
    yr <- y
  }
  yr <- yr[order(yr[,1]),]
  colnames(xr1) <- colnames(yr) <- colnames(x)

  # cut down to elements of xr in yr
  xr <- xr1[xr1$taxon %in% yr$taxon,]
  rownames(xr) <- NULL
  # dataframe to store output
  z <- data.frame(xr$taxon, "000", NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(z) <- c("taxon", "status", "n_1R1", "n_0R1", "n_1R0", "n_0R0", "n_00R", "n_R00", "fad_diff", "lad_diff")
  flag <- rep("000", times = nrow(x))
  fad_diff <- flag
  fad_diff[] <- NA
  lad_diff <- fad_diff

  # convert to data.table for fast indexing
  x$rownum <- 1:nrow(x)
  x <- data.table::as.data.table(x)
  data.table::setkeyv(x, "taxon")

  # for each unique genus name in x
  for(i in 1:nrow(xr)) {

    # get all occurrences of the taxon
    occs <- x[.(xr$taxon[i]),]
    # get the matching range in the checking dataset
    mt <- which(yr$taxon == xr$taxon[i])
    # max and min ages
    upr <- max(occs$max)
    lwr <- min(occs$min)

    # if the taxon is not present in the checking dataframe, assign skip
    if(length(mt) == 0) {
      1 + 1
    } else {

      z[i,3:8] <- 0
      # if all ranges are valid
      if(all(occs$max <= yr$max[mt]) & all(occs$min >= yr$min[mt])) {
        flag[occs$rownum] <- "R1R"
        z$n_1R1[i] <- nrow(occs)
        z$status[i] <- "R1R"
        z$fad_diff[i] <- z$lad_diff[i] <- NA

      } else {

        # check for any valid occurrences
        v <- occs[which(occs$max <= yr$max[mt] &
                          occs$min >= yr$min[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "R1R"
          z$n_1R1[i] <- nrow(v)
          occs <- occs[-(which(occs$max <= yr$max[mt] &
                                 occs$min >= yr$min[mt])),]
        }

        # check for occurrences which exceed FAD and LAD
        v <- occs[which(occs$min < yr$min[mt] &
                          occs$max > yr$max[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "0R0"
          fad_diff[v$rownum] <- v$max - yr$max[mt]
          lad_diff[v$rownum] <- yr$max[mt] - v$min
          z$n_0R0[i] <- nrow(v)
          z$status[i] <- "0R0"
        }

        # check for occurrences fully outside of the FAD or LAD
        v <- occs[which(occs$max > yr$max[mt] &
                          occs$min > yr$max[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "00R"
          fad_diff[v$rownum] <- v$max - yr$max[mt]
          z$n_00R[i] <- nrow(v)
          if(z$status[i] == "000") {z$status[i] <- "00R"}
        }
        v <- occs[which(occs$min < yr$min[mt] &
                          occs$max < yr$min[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "R00"
          lad_diff[v$rownum] <- yr$max[mt] - v$min
          z$n_R00[i] <- nrow(v)
          if(z$status[i] == "000") {z$status[i] <- "R00"}
          if(z$status[i] == "00R") {z$status[i] <- "0R0"}
        }

        # check for occurrences which exceed the FAD only
        v <- occs[which(occs$min <= yr$max[mt] & occs$min >= yr$min[mt] &
                          occs$max > yr$max[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "0R1"
          fad_diff[v$rownum] <- v$max - yr$max[mt]
          z$n_0R1[i] <- nrow(v)
          if(z$status[i] == "000") {z$status[i] <- "0R1"}
          if(z$status[i] == "R00") {z$status[i] <- "0R0"}
        }
        # check for occurrences which exceed the LAD only
        v <- occs[which(occs$min < yr$min[mt] &
                          occs$max >= yr$min[mt] & occs$max <= yr$max[mt]),]
        if(nrow(v) != 0) {
          flag[v$rownum] <- "1R0"
          lad_diff[v$rownum] <- yr$max[mt] - v$min
          z$n_1R0[i] <- nrow(v)
          if(z$status[i] == "000") {z$status[i] <- "1R0"}
          if(z$status[i] %in% c("00R", "01R")) {z$status[i] <- "0R0"}
        }

        # range difference checks
        if(upr > yr$max[mt]) {
          z$fad_diff[i] <- upr - yr$max[mt]
        } else {
          z$fad_diff[i] <- NA
        }
        if(lwr < yr$min[mt]) {
          z$lad_diff[i] <- yr$min[mt] - lwr
        } else {
          z$lad_diff[i] <- NA
        }

      }
    }
    # notify R
    if(verbose) {
      if(i != 1) {cat(paste0("\r"))}
      cat(paste0("Taxon ", i, "/", nrow(xr), " checked"))
      if(i == nrow(xr)) {cat("\n")}
    }
  }

  ftails <- do.call(cbind, lapply(flag.diff, function(x) {
    tail_logical <- as.numeric(fad_diff >= x)
    tail_logical[is.na(tail_logical)] <- 0
    tail_logical
  }))
  ltails <- do.call(cbind, lapply(flag.diff, function(x) {
    tail_logical <- as.numeric(lad_diff >= x)
    tail_logical[is.na(tail_logical)] <- 0
    tail_logical
  }))
  per_occ <- cbind.data.frame(flag, round(as.numeric(fad_diff), digits = 2),
                              round(as.numeric(lad_diff), digits = 2), ftails, ltails)
  colnames(per_occ) <- c("status", "fad_diff", "lad_diff", paste0("fad_flag.", flag.diff), paste0("lad_flag.", flag.diff))

  # return
  if(verbose) {message("See $occurrence in output for the error statuses of individual occurrences")}
  out <- list()
  out[[1]] <- z
  out[[2]] <- per_occ
  names(out) <- c("taxon", "occurrence")
  return(out)
}
