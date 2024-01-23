#' threshold_peaks
#'
#' Function to detect if two peaks in a density spectrum can
#' be considered separate based on a user supplied threshold.
#' Creates a sequence of divisions from the troughs immediately
#' preceding any significant peaks, then bins occurrences for
#' a given taxon name by those divisions.
#' @param x A list of significant peaks as returned by `find_peaks`
#' @param y An occurrence dataset with taxon names corresponding
#' to the list names of x
#' @param ycols A character vector denoting, in order, the taxon, FAD
#' and LAD columns in y
#' @param thresh The threshold distance between peaks above which
#' they will be considered distinct - given in Ma
#' @param verbose A logical determining if function progress
#' should be reported
#' @import data.table

threshold_peaks <- function(x, y, ycols = c("genus", "max_ma", "min_ma"), thresh = 15, verbose = TRUE) {

  # check args
  if(!is.list(x)) {
    stop("x should be a list of named numeric vectors")
  }
  if(!all(is.numeric(unlist(lapply(x, class))))) {
    stop("x should be a list of named numeric vectors")
  }
  if (!is.data.frame(y)) {
    stop("y must both be dataframes")
  }
  if (ncol(y) < 3) {
    stop("Supplied dataframes do not contain enough data (at least name, FAD and LAD needed in both")
  }
  if (length(ycols) != 3) {
    stop("ycols must contain three elements")
  }
  if (!all(ycols %in% colnames(y))) {
    stop("One or more elements of ycols are not column names in y")
  }
  if (!is.numeric(y[, ycols[2]]) | !is.numeric(y[, ycols[3]])) {
    stop("Elements 2 and 3 of ycols must refer to numeric columns in y")
  }
  y <- y[,ycols]
  if(length(thresh) != 1 && length(thresh) != length(x)) {
    stop("thresh must be either a single positive numeric which will be applied to all
         elements of x, or the same length of x")
  }
  if(any(thresh < 1)) {
    stop("All elements of thresh must be positive integers")
  }
  thresh <- round(thresh)
  if(length(thresh) == 1) {
    thresh <- rep(thresh, times = length(x))
  }

  # set up additional elements
  . <- NULL
  y2 <- cbind.data.frame(y, y[,2] - ((y[,2] - y[,3]) / 2), 1:nrow(y))
  colnames(y2) <- c("taxon", "max_ma", "min_ma", "mpt", "rownum")
  y2 <- data.table::as.data.table(y2)
  data.table::setkeyv(y2, "taxon")

  # set up storage
  cat_out <- rep(NA, nrow(y2))
  sig_peaks <- list()

  # scan through
  for(i in 1:length(x)) {

    # get occs
    occs <- y2[.(names(x)[i]),]

    # get peak trough sequence
    foo <- x[[i]]

    # if there is only one peak or none
    if(length(foo) < 3) {
      sig_pk <- foo
      divs <- c(ceiling(max(occs$max_ma)), floor(min(occs$min_ma)))
      cat_out[occs$rownum] <- paste("1/1")

    } else {

      # get peaks
      pk <- foo[seq(from = 1, to = length(foo), by = 2)]
      # difference between peaks
      pdiff <- abs(diff(as.numeric(names(pk))))
      # apply threshold
      sdiff <- pdiff >= thresh[i]

      # if all the peaks collapse together
      if(sum(sdiff) == 0) {
        divs <- c(ceiling(max(occs$max_ma)), floor(min(occs$min_ma)))
        cat_out[occs$rownum] <- paste("1/1")

        # otherwise get the boundaries
      } else {
        # get initial threshold, setting the first peak as true
        pseq <- c(1, as.numeric(sdiff))
        # get the cumulative sum
        pseq2 <- cumsum(pseq)
        names(pseq2) <- names(pk)
        # get the unique values (i.e. where cumsum has detected TRUE)
        pseq3 <- unique(pseq2)
        # set the last value to the last peak so that the bounds are always included
        pseq3 <- c(pk[pseq3[-length(pseq3)]], pk[length(pk)])
        # as unique is selected forwards, the troughs are selected backwards so that
        # trough immediately preceding the next unique peak is taken
        tseq <- match(names(pseq3), names(foo))
        tseq <- foo[tseq[2:length(tseq)] - 1]
        sig_pk <- as.vector(na.omit(as.vector(rbind(pseq3, c(tseq, NA)))))
        names(sig_pk) <- as.vector(na.omit(as.vector(rbind(names(pseq3), c(names(tseq), NA)))))
        # revision - return troughs as these are the true boundaries to be used
        divs <- c(ceiling(max(occs$max_ma)), as.numeric(names(tseq)), floor(min(occs$min_ma)))
        # slice
        cuts <- cut(occs$mpt, breaks = divs, labels = FALSE)
        cat_out[occs$rownum] <- paste0(cuts, "/", max(cuts))
      }
    }
    sig_peaks[[i]] <- divs

    # notify R
    if(verbose) {
      if(i != 1) {cat(paste0("\r"))}
      cat(paste0("Taxon ", i, "/", length(x), " checked"))
    }
  }
  names(sig_peaks) <- names(x)

  # compare lengths
  c1 <- unlist(lapply(x, function(z) {
    floor((length(z) + 1) / 2)
  }))
  c2 <- unlist(lapply(sig_peaks, function(z) {
    length(z) - 1
  }))
  c3 <- cbind(c1, c2)
  colnames(c3) <- c("pre", "post")
  rownames(c3) <- names(x)

  # output
  out <- list(sig_peaks, cat_out, c3)
  names(out) <- c("divisions", "groups", "counts")
  return(out)
}
