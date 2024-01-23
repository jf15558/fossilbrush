#' threshold_ranges
#'
#' Function to detect if two peaks in a density spectrum can
#' be considered separate based on a user supplied threshold.
#' Creates a sequence of divisions from the troughs immediately
#' preceding any significant peaks, then bins occurrences for
#' a given taxon name by those divisions.
#' @param x An occurrence dataset containing taxon names,
#' maximum ages and minimum ages
#' @param rank The column name in x containing the taxon names
#' @param srt A column name in x denoting the occurrence maximum ages
#' @param end A column name in x denoting the occurrence minumum ages
#' @param method The method for quantifying occurrence density: one
#' histogram or kernel. Kernel is the recommended default. As called
#' be @seealso densify
#' @param step A positive integer specifying the time window size
#' for density calculation. As called by @seealso densify
#' @param density A positive numeric specifying the step size for
#' densifying records. This should ideally be smaller than step.
#' As called by @seealso densify
#' @param use_sd A logical determining whether to use peaks detected
#' as significant using the mean + standard deviation of its neighbourhood.
#' If FALSE, then the peaks need only be greater than the neighbourhood
#' mean to be significant. Thus, use_sd is more conservative, but less
#' prone to noise. As called by @seealso find_peaks
#' @param win A positive integer specifying the neighborhood window length
#' on either side of a peak durign significance testing (i.e. win 5 will
#' give a total window of 11: -5 indices + peak index + 5 indices). As
#' called by @seealso find_peaks
#' @param thresh The threshold distance between peaks above which
#' they will be considered distinct - given in Ma
#' @param ... additional arguments passed to @seealso density
#' @param report A logical determining if the analytical outputs of the
#' function be returned to the user, as well as the revised taxon names, TRUE
#' by default
#' @param verbose A logical determining if function progress
#' should be reported
#' @return If report = TRUE (the default), a list of five elements. $data
#' gives the thresholded (and potentially subdivided) taxon names. $matrix
#' is the taxon-wise matrix of occurrence densities. $peaks is a list
#' containing three lists of peaks (all peaks, significant by mean + sd,
#' significant by sd only) for each taxon and a dataframe of peak counts
#' between the three treatments. $comparison
#' @import data.table
#' @import pbapply
#' @import data.table
#' @importFrom methods as slot
#' @importFrom stats na.omit supsmu sd
#' @importClassesFrom Matrix sparseVector
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # interpeak thresholding
#' itp <- threshold_ranges(brachios, win = 8, thresh = 10,
#'                         rank = "genus", srt = "max_ma", end = "min_ma")

threshold_ranges <- function(x, rank = "genus", srt = "max_ma", end = "min_ma", method = "kernel",
                            step = 1, density = 0.1, use_sd = TRUE, win = 5, thresh = 5, ..., report = TRUE, verbose = TRUE) {

  # data supply args
  if(!exists("x")) {
    stop("Please supply x as a dataframe of occurrences with taxonomic names, and maximum and minimum ages ages")
  }
  if(!is.data.frame(x)) {
    stop("Please supply x as a dataframe of occurrences with taxonomic names, and maximum and minimum ages ages")
  }
  if(ncol(x) < 3) {
    stop("Please supply x as a dataframe of occurrences with taxonomic names, and maximum and minimum ages ages")
  }
  if(!is.character(rank) | length(rank) != 1) {
    stop("rank should be a character vector of length one specifing the taxonomic names column in x")
  }
  if(!is.character(srt) | length(srt) != 1) {
    stop("srt should be a character vector of length one specifing the maximum ages column in x")
  }
  if(!is.character(end) | length(end) != 1) {
    stop("end should be a character vector of length one specifing the maximum ages column in x")
  }
  if(!all(c(rank, srt, end) %in% colnames(x))) {
    stop("One or more of rank, srt or end are not colnames in data")
  }
  if(!is.character(x[,rank])) {
    stop("Column rank must be character")
  }
  if(!is.numeric(x[,srt]) | !is.numeric(x[,end])) {
    stop("Columns srt and end must be numeric")
  }
  if(any(c(is.na(x[,srt]), is.na(x[,end])))) {
    stop("Missing values are present in columns srt and/or end")
  }
  if(any(x[,srt] < x[,end])) {
    stop("One or more maximum ages are smaller than their corresponding minimum ages")
  }

  # check analytical args (aside from thresh)
  if(!is.character(method) | length(method) != 1) {
    stop("method should be a character of length 1: one of histogram or kernel")
  }
  if(!all(method %in% c("histogram", "kernel"))) {
    stop("Method must be one or both of histogram or kernel")
  }
  if(!is.numeric(step) | length(step) != 1) {
    stop("step should be a numeric of length 1 ")
  }
  if(step < 1 | step %% 1 != 0) {
    stop("Step must be a single positive integer")
  }
  if(!is.numeric(density) | length(density) != 1) {
    stop("density should be a numeric of length 1 ")
  }
  if(step < 0) {
    stop("density must be positive")
  }
  if(density >= step) {
    warning("Density should ideally be smaller than step")
  }
  if(!is.numeric(win) | length(win) != 1) {
    stop("win should be a numeric of length 1 ")
  }
  if(win < 1 | win %% 1 != 0) {
    stop("win must be a single positive integer")
  }
  if(!is.logical(verbose) | length(verbose) != 1) {
    stop("verbose should be a logical of length 1")
  }

  # copy then trim data
  x1 <- x
  x <- cbind.data.frame(x[,rank], x[,srt], x[,end], x[,srt] - ((x[,srt] - x[,end]) / 2), 1:nrow(x))
  colnames(x) <- c(rank, srt, end, "mpt", "rownum")
  x <- x[complete.cases(x),]
  # add a very small constant to deal with zero range ages
  x[,srt] <- x[,srt] + density
  # set the upper and lower range for the vectors
  start_v <- ceiling(max(x[,srt]))
  end_v <- floor(min(x[,end]))
  x <- data.table::as.data.table(x)
  data.table::setkeyv(x, rank)

  # check thresh arg
  tax <- unique(x[[rank]])
  if(is.atomic(thresh)) {
    if(!is.numeric(thresh) | length(thresh) != 1) {
      stop("If not supplying a dataframe, thresh must be a single positive numeric. This value will be used for all taxa")
    }
    if(thresh <= 0) {
      stop("If not supplying a dataframe, thresh must be a single positive numeric. This value will be used for all taxa")
    }
    tax <- tax[order(tax)]
    thresh <- cbind.data.frame(tax, rep(thresh, length(tax)))

  } else {
    if(!is.data.frame(thresh)) {
      stop("To supply taxon-specific thresholds, thresh should be a dataframe with a column of taxon names and a column of taxon thresholds")
    }
    if(ncol(thresh) != 2) {
      stop("To supply taxon-specific thresholds, thresh should be a dataframe with a column of taxon names and a column of taxon thresholds")
    }
    if(nrow(thresh) != length(tax)) {
      stop("To supply taxon-specific thresholds, thresh should be a dataframe with as many rows as unique taxa (i.e unique elements in x[,ranks])")
    }
    if(!is.character(thresh[,1])) {
      stop("The first column of thresh should be a character vector of taxon names")
    }
    if(!all(thresh[,1] %in% tax)) {
      stop("Not all taxon names in x have a corresponding entry in thresh")
    }
    if(!is.numeric(thresh[,2])) {
      stop("The second column of thresh should be a numeric vector of taxon-specific thresholds")
    }
    if(any(thresh <= 0)) {
      stop("Supplied thresholds must all be positive")
    }
    thresh <- thresh[order(thresh[,1]),]
  }
  colnames(thresh) <- c("tax", "thresh")

  # global variable workaround for data.table
  . <- slot <- NULL

  # internally define sv_cbind function as is very small
  sv_cbind <- function (...) {
    input <- lapply(list(...), as, "dsparseVector")
    thelength <- unique(sapply(input,length))
    stopifnot(length(thelength) == 1)
    return(Matrix::sparseMatrix(
      x = unlist(lapply(input, slot, "x")),
      i = unlist(lapply(input, slot, "i")),
      p = c(0, cumsum(sapply(input, function(x) {length(x@x)}))),
      dims = c(thelength, length(input))
    ))
  }

  # for each taxon
  tp <- pbsapply(1:length(tax), simplify = FALSE, function(y) {

    ######## DENSITY ########

    # get all occurrences of the taxon (uses data.table for fast indexing)
    upr_occ_m <- x[.(tax[y]),]
    upr_occ <- upr_occ_m[,2:3]
    # sequence from each FAD-LAD pair by density
    upr <- as.vector(unlist(apply(upr_occ, 1, function(z) {seq(from = z[1], to = z[2], by = -density)})))
    # set the bins for the density calculation
    start_range <- ceiling((max(upr)))
    end_range <- floor((min(upr)))
    bins <- seq(from = end_range, to = start_range, by = step)

    # set additional vectors to pad around the ranges
    if(start_range != start_v) {
      pre_base <- rep(0, length(seq(from = start_v, to = start_range + step, by = -step)))
    } else {
      pre_base <- NULL
    }
    if(end_range != end_v) {
      post_top <- rep(0, length(seq(from = end_range - step, to = end_v, by = -step)))
    } else {
      post_top <- NULL
    }

    # get the density of the 'densified' records per time slice
    upr_h <- NULL
    if("histogram" %in% method | "kernel" %in% method) {
      # do the histogram and divide by the densification factor (i.e. extra points created per bin), ceiling to be safe
      upr_h <- ceiling(c(pre_base, rev(hist(upr, breaks = bins, plot = FALSE)$counts), post_top) / (step / density))
    }
    foo <- NULL
    if("kernel" %in% method) {

      # if only 1 bin, put all density in that bin
      if(length(bins) < 3) {
        foo <- c(pre_base, 1, post_top)
      } else {
        # otherwise do density across bins
        den <- rev(density(upr, n = (length(bins) - 1), from = end_range, to = start_range, ...)$y)
        # if statement to catch very rare 0 density estimates (unclear why, but only 6 cases out of 66000)
        if(sum(den) == 0) {
          den <- rep(1 / length(bins - 1), times = length(bins) - 1)
        }
        foo <- c(pre_base, den, post_top)
        foo[which(upr_h == 0)] <- 0
      }
    } else {
      foo <- upr_h
    }
    names(foo) <- seq(from = start_v, to = end_v + step, by = -step)
    # bind and store as sparse vector
    upr2 <- as(foo, "sparseVector")

    ######## PEAKS ########

    # get peaks
    foo2 <- foo[c(TRUE, (diff(foo) != 0))]
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
      vec <- vec2 <- vector()
      for(i in 1:length(pk)) {

        # get the numbers within the window, padding with zeros if window extends beyond data time extent
        winrng <- c((match(names(pk[i]), names(foo)) - win), (match(names(pk[i]), names(foo)) + win))
        winrng2 <- winrng
        if(winrng[1] < 1) {winrng[1] <- 1}
        if(winrng[2] > length(foo)) {winrng[2] <- length(foo)}
        winrng3 <- abs(winrng - winrng2)
        addz1 <- rep(0, winrng3[1])
        addz2 <- rep(0, winrng3[2])
        vals <- c(addz1, foo[winrng[1]:winrng[2]], addz2)

        #vals <- foo[(match(names(pk[i]), names(foo)) - win):(match(names(pk[i]), names(foo)) + win)]
        # extra smoothing
        vals <- as.vector(na.omit(vals))
        vals2 <- supsmu(1:length(vals), vals, span = 0.2)$y

        if(max(vals2[win:(win + 2)]) > ((mean(na.omit(vals2))) + ((sd(na.omit(vals2)))))) {
          vec[i] <- TRUE
        } else {
          vec[i] <- FALSE
        }
        if(max(vals2[win:(win + 2)]) > (mean(na.omit(vals2)))) {
          vec2[i] <- TRUE
        } else {
          vec2[i] <- FALSE
        }
      }
      # if all peaks are insigificant, take the maximum as the single valid peak (bit of fudge), but this is v unlikely to be needed
      if(sum(vec) == 0) {vec <- which.max(pk)}
      if(sum(vec2) == 0) {vec2 <- which.max(pk)}

      # select the thresholded peaks
      pt_ms <- as.vector(pt[,vec])
      # remove the tail value as it is the unimportant trough (this works for any peak number)
      pt_ms <- pt_ms[-length(pt_ms)]
      # retrieve names
      ptn_ms <- as.vector(ptn[,vec])
      names(pt_ms) <- ptn_ms[-length(ptn_ms)]

      # select the thresholded peaks
      pt_m <- as.vector(pt[,vec2])
      # remove the tail value as it is the unimportant trough (this works for any peak number)
      pt_m <- pt_m[-length(pt_m)]
      # retrieve names
      ptn_m <- as.vector(ptn[,vec2])
      names(pt_m) <- ptn_m[-length(ptn_m)]
    }
    out_apply <- list(pt_all, pt_ms, pt_m)

    ######## THRESHOLD ########

    # get peak trough sequence
    if(use_sd) {
      use_peak <- pt_ms
    } else {
      use_peak <- pt_m
    }

    # if there is only one peak or none
    if(length(use_peak) < 3) {
      sig_pk <- use_peak
      divs <- c(ceiling(max(upr_occ_m[[srt]])), floor(min(upr_occ_m[[end]])))
      cat_out <- cbind(upr_occ_m[["rownum"]], "1/1")

    } else {

      # get peaks
      pk2 <- use_peak[seq(from = 1, to = length(use_peak), by = 2)]
      # difference between peaks
      pdiff <- abs(diff(as.numeric(names(pk2))))
      # apply threshold (this covers peak 2 to peak n, as diff will not include peak 1)
      sdiff <- pdiff >= thresh[y,2]

      # if all the peaks collapse together
      if(sum(sdiff) == 0) {
        divs <- c(ceiling(max(upr_occ_m[[srt]])), floor(min(upr_occ_m[[end]])))
        cat_out <- cbind(upr_occ_m[["rownum"]], "1/1")

        # otherwise get the boundaries
      } else {

        # sdiff applies to peak 2 onwards, so [-1] to exclude first peak from the logical indexing
        pseq <- names(pk2)[-1][sdiff]
        # match the significant peak names back into the peak trough sequence
        tseq <- match(pseq, names(use_peak))
        # matched peak indices - 1 to get their preceding troughs
        tseq <- use_peak[tseq - 1]
        # add in the ends of the total name range (i.e the trough for the first peak, automatically included as the first peak must be valid, and the last trough as is the range terminus)
        divs <- c(ceiling(max(upr_occ_m[[srt]])), as.numeric(names(tseq)), floor(min(upr_occ_m[[end]])))
        # slice
        cuts <- cut(upr_occ_m[["mpt"]], breaks = divs, labels = FALSE)
        cat_out <- cbind(upr_occ_m[["rownum"]], paste0(cuts, "/", max(cuts)))
      }
    }
    split_data <- list(divs, cat_out)

    # set up pb output
    p_out <- list()
    # density matrix column
    p_out[[1]] <- upr2
    # peak index vectors
    p_out[[2]] <- out_apply
    # name split data
    p_out[[3]] <- split_data

    # return from pb
    return_pb <- p_out
  })

  # get density matrix columns
  dm <- lapply(tp, function(x) {x[[1]]})
  # bind into matrix
  dm <- do.call(sv_cbind, dm)
  # name dimensions
  dimnames(dm) <- list(seq(from = start_v, to = end_v + step, by = -step), tax)

  # get peak data
  pks <- lapply(tp, function(x) {x[[2]]})
  # extract out the all, standard deviation and no standard deviation series
  pks_all <- lapply(pks, function(z) {z[[1]]})
  pks_ms <- lapply(pks, function(z) {z[[2]]})
  pks_m <- lapply(pks, function(z) {z[[3]]})
  # count peaks in each series
  c1 <- unlist(lapply(pks_all, function(z) {floor((length(z) + 1) / 2)}))
  # count ms peaks (troughs will automatically be n(peaks) - 1)
  c2 <- unlist(lapply(pks_ms, function(z) {floor((length(z) + 1) / 2)}))
  # count m peaks (troughs will automatically be n(peaks) - 1)
  c3 <- unlist(lapply(pks_m, function(z) {floor((length(z) + 1) / 2)}))
  # make dataframe and name elements
  counts <- cbind.data.frame(c1, c2, c3)
  colnames(counts) <- c("all", "mean_sd", "mean")
  rownames(counts) <- tax
  # finalise
  names(pks_all) <- names(pks_ms) <- names(pks_m) <- tax
  pks_out <- list(pks_all, pks_ms, pks_m, counts)
  names(pks_out) <- c("profile_all", "profile_ms", "profile_m", "counts")

  # split data
  splits <- lapply(tp, function(x) {x[[3]]})
  # diversity
  c1 <- unlist(lapply(splits, function(z) {
    floor((length(z[[1]]) + 1) / 2)
  }))
  c2 <- unlist(lapply(splits, function(z) {
    length(z[[1]]) - 1
  }))
  c3 <- cbind(c1, c2)
  colnames(c3) <- c("pre", "post")
  rownames(c3) <- tax
  # names
  div <- lapply(splits, function(x) {x[[2]]})
  div <- do.call(rbind.data.frame, div)
  gen2 <- x1[,rank]
  gen2[as.numeric(div[,1])] <- paste(gen2[as.numeric(div[,1])], div[,2], sep = "_")

  # output
  out <- list(gen2, dm, pks_out, c3)
  names(out) <- c("data", "matrix", "peaks", "comparison")
  if(!report) {
    out <- gen2
  }
  return(out)
}
