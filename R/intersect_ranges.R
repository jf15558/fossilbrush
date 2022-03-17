#' intersect_ranges
#'
#' Function to find the maximum intersection between a set
#' of numeric ranges, in this case first and last appearence
#' datums on taxonomic ranges.
#' @param x A numeric data.frame or matrix of ranges. If
#' just two columns are supplied, the first column is assumed
#' to be the srt column
#' @param srt If x contains more than two columns, srt is the
#' name of the range base column - the FAD
#' @param end If x contains more than two columns, end is the
#' name of the range top column - the LAD
#' @param verbose A logical indicating whether the function
#' should report progress to the console
#' @return A matrix with three columns, indicating the intersection
#' (FAD and LAD) and the number of ranges that intersection
#' encompasses
#' @export

intersect_ranges <- function(x, srt = NULL, end = NULL, verbose = TRUE) {

  # check x
  if(!exists("x")) {
    stop("Please supply x as a dataframe or matrix of ranges")
  }
  if(!class(x)[1] %in% c("data.frame", "matrix")) {
    stop("Please supply x as a dataframe or matrix of ranges")
  }
  if(ncol(x) < 2) {
    stop("x must minimally contain a srt column and a lad column")
  }
  if(ncol(x) > 2) {
    if(is.null(srt) | is.null(end))
      stop("If x contains more than two columns, srt and end must be specified")
    x <- x[, c(srt, end)]
  }
  if(!is.numeric(x[,1]) | !is.numeric(x[,2])) {
    stop("x must be numeric")
  }
  if(any(is.na(x))) {
    warning("x contains missing values - these rows will be dropped")
    x <- x[complete.cases(x),]
  }
  if(nrow(x) < 1) {
    stop("x does not contain any rows (this may be due to internal removal of incomplete rows)")
  }
  if(any(x[,srt] < x[,end])) {
    stop("One or more maximum ages in x are smaller than their corresponding minimum ages")
  }

  # rank everything
  # in case the same date is given, they will have the same rank.
  # this will create additional 'slots in the logical matrix,
  # which should not be an issue.
  ranked <- rank(x, ties.method = "min")

  # copy the dimensions
  dim(ranked) <- dim(x)

  # compare the ranks
  where <- rep(0, max(ranked) * 2 - 1)
  for(i in 1:nrow(x)){
    # index to increment
    ind <- seq(ranked[i,1] * 2 - 1, ranked[i,2] * 2 - 1)
    # increment in this interval
    where[ind] <- where[ind] + 1
  }

  # no overlap
  if(all(where == 1)){
    if(verbose)	message("No overlap found")
    results <- cbind(1, 1, 1)
    colnames(results) <- c("lb", "ub", "N")
  } else {

    # find the indices!
    # what is the maximum overlap?
    maxVals <- max(where)

    # print to the terminal
    if(verbose) {
      # report the number of overlapping ranges
      if(maxVals == nrow(x)){
        message("All ranges overlap.")
      } else {
        message(paste(maxVals, "out of the", nrow(x), "ranges overlap."))
      }
    }

    # are there multiple solutions???
    streaks <- which(maxVals == where)

    # are there multiple solutions???
    # shift the index of the difference by 1
    diffs <- c(0, diff(streaks))

    # changepoints between solutions (index of diffs)
    # the number of solutions is the length
    changes <- c(1, which(diffs > 1))

    # final result - one row for every solution!
    results <- matrix(ncol = 3, nrow = length(changes))

    # for every solution
    for(i in 1:length(changes)){
      # normal case
      if(i != length(changes)){
        # one solutions
        solStreaks <- streaks[changes[i]:(changes[i + 1] - 1)]
        # in case there are ties, there will be overlaps that leads to repetitions
        # last iteration
      } else {
        solStreaks <- streaks[changes[i]:length(streaks)]
      }

      # find the values that correspond to the indices
      fadInLogical <- min(solStreaks)
      ladInLogical <- max(solStreaks)

      # has to be transformed
      # do an assertion
      # the overlap intervals should range from date to date.
      # neither of these should be even col. indices in the logical matrix,
      # or things are f-ed up.
      if(fadInLogical%%2 == 0) stop("FAD in logical matrix should not be even!")
      if(ladInLogical%%2 == 0) stop("FAD in logical matrix should not be even!")

      # the rank of starting interval
      startDateRank <- (fadInLogical + 1) / 2
      endDateRank <- (ladInLogical + 1) / 2

      # corresponding dates
      first <- x[ranked == startDateRank]
      last <- x[ranked == endDateRank]

      # take the first result in the case of potential ties that rank() outputs.
      results[i, ] <- c(first[1],last[1], maxVals[1])
    }
    # copy the column names of the original entry
    colnames(results) <- c("lb", "ub", "N")
    # indicate how many overlaps there are
    rownames(results) <- paste0(maxVals, "_", 1:nrow(results))
  }
  # return object
  return(as.data.frame(results))
}
