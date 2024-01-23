#' find_peaks
#'
#' Function to scan, column-wise, a matrix of per-taxon observation
#' density time series. This can be applied to either the histogram
#' or the kernel density output of `densify`, but the latter is
#' recommended. Peaks are detected as local maxima, then smoothed
#' within a local window and tested to distinguish if they are
#' noise or significant. Strict threshold is that the peak is
#' greater than the mean + sd of the window
#' @param x A matrix as outputed by `densify`
#' @param win A positive integer specifying the window length on
#' either side of a peak (i.e. win 5 will give a total window of 11 -
#' -5 indices + peak index + 5 indices)
#' @param verbose A logical determining if function progress
#' should be reported
#' @return A list of four, the first three positions containing lists
#' of the peak indices for each taxon, under raw, mean + sd and mean
#' detection regimes. The fourth item is a dataframe of counts of peaks
#' per taxon, 1 row per taxon, 1 column per detection regime
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom stats na.omit supsmu sd
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # get density matrix
#' dens <- densify(brachios)
#' # run function, using kernel density matrix
#' pk <- find_peaks(dens$kdensity)

find_peaks <- function(x, win = 5, verbose = TRUE) {

  if(!exists("x")) {
    stop("Please supply x as a matrix as outputted by densify")
  }
  if(!class(x)[1] %in% c("matrix", "dgCMatrix", "Matrix")) {
    stop("Please supply x as a matrix as outputted by densify")
  }
  if(length(win) != 1) {
    stop("win must be a single positive integer")
  }
  if(win < 1) {
    stop("win must be a single positive integer")
  }
  if(!verbose) {
    baseopt <- getOption("pboptions")
    opb <- pboptions(type = "none")
  }
  win <- round(win)

  out1 <- pbsapply(1:ncol(x), simplify = FALSE, function(y) {

    foo <- x[,y]
    names(foo) <- rownames(x)
    foo2 <- foo[c(TRUE, (diff(foo) != 0))]

    # get peaks
    a1 <- c(0, foo2, 0)
    a2 <- c(foo2, 0, 0)
    a3 <- c(0, 0, foo2)
    pk <- foo2[match(names(which((a1 > a2 & a1 > a3)[2:(length(a1) - 1)])), names(foo2))]
    # get troughs
    tr <- foo2[match(names(which((a1 < a2 & a1 < a3)[2:(length(a1) - 1)])), names(foo2))]
    # bind together into peak-trough sequence (na used to pad to even lengths)
    pt <- rbind(pk, c(tr, NA))
    ptn <- rbind(names(pk), c(names(tr), NA))
    pt_all <- as.vector(na.omit(as.vector(pt)))
    names(pt_all) <- as.vector(na.omit(as.vector(ptn)))

    # assess peak validity if multiple peaks
    if(length(pt) == 2) {
      pt_final <- pt[1]
      names(pt_final) <- colnames(pt)
      pt_ms <- pt_m <- pt_final

    } else {

      # assess peak validity within its neighbourhood in a given window
      # either side of the peak
      vec <- vector()
      for(i in 1:length(pk)) {
        vals <- foo[(match(names(pk[i]), names(foo)) - win):(match(names(pk[i]), names(foo)) + win)]
        fv <- vals[win + 1]
        vals <- as.vector(na.omit(vals))
        # extra smoothing
        vals2 <- supsmu(1:length(vals), vals)$y
        if(vals2[win + 1] > ((mean(na.omit(vals2))) + ((sd(na.omit(vals2)))))) {
          vec[i] <- TRUE
        } else {
          vec[i] <- FALSE
        }
      }
      # select the thresholded peaks
      pt_ms <- as.vector(pt[,vec])
      # remove the tail value as it is the unimportant trough (this works for any peak number)
      pt_ms <- pt_ms[-length(pt_ms)]
      # retrieve names
      ptn_ms <- as.vector(ptn[,vec])
      names(pt_ms) <- ptn_ms[-length(ptn_ms)]

      # redo without the standard deviation adjustment
      for(i in 1:length(pk)) {
        vals <- foo[(match(names(pk[i]), names(foo)) - win):(match(names(pk[i]), names(foo)) + win)]
        fv <- vals[win + 1]
        vals <- as.vector(na.omit(vals))
        # extra smoothing
        vals2 <- supsmu(1:length(vals), vals)$y
        if(vals2[win + 1] > (mean(na.omit(vals2)))) {
          vec[i] <- TRUE
        } else {
          vec[i] <- FALSE
        }
      }
      # select the thresholded peaks
      pt_m <- as.vector(pt[,vec])
      # remove the tail value as it is the unimportant trough (this works for any peak number)
      pt_m <- pt_m[-length(pt_m)]
      # retrieve names
      ptn_m <- as.vector(ptn[,vec])
      names(pt_m) <- ptn_m[-length(ptn_m)]
    }
    out_apply <- list(pt_all, pt_ms, pt_m)
  })
  out_all <- lapply(out1, function(z) {
    z[[1]]
  })
  out_ms <- lapply(out1, function(z) {
    z[[2]]
  })
  out_m <- lapply(out1, function(z) {
    z[[3]]
  })

  # count all peaks
  c1 <- unlist(lapply(out1, function(z) {
    floor((length(z[[1]]) + 1) / 2)
  }))
  # count ms peaks (troughs will automatically be n(peaks) - 1)
  c2 <- unlist(lapply(out1, function(z) {
    floor((length(z[[2]]) + 1) / 2)
  }))
  # count m peaks (troughs will automatically be n(peaks) - 1)
  c3 <- unlist(lapply(out1, function(z) {
    floor((length(z[[3]]) + 1) / 2)
  }))
  counts <- cbind.data.frame(c1, c2, c3)
  colnames(counts) <- c("all", "mean_sd", "mean")
  rownames(counts) <- colnames(x)

  # format output
  names(out_all) <- names(out_ms) <- names(out_m) <- colnames(x)
  out <- list(out_all, out_ms, out_m, counts)
  names(out) <- c("profile_all", "profile_ms", "profile_m", "counts")
  if(!verbose) {opt <- pboptions(baseopt)}
  return(out)
}
