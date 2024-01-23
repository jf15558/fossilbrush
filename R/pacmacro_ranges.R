#' pacmacro_ranges
#'
#' Function to apply a modification of Pacman trimming to
#' macrofossil data. The function generates a densified
#' occurrence record using the same methods as `densify`
#' then trim the upper and lower ranges by a user-defined
#' percentage. The full and trimmed ranges are then
#' compared against each other to test if the FAD and the
#' LAD for a taxon form a long tail in its distribution.
#' Multiple tail thresholds can be supplied, but all test
#' to see if the sum of the FAD and LAD which exceeds the
#' trimmed range constitute the threshold proportion of the
#' total range for than taxon, e.g. does the FAD and the
#' LAD outside of the trimmed range comprise a quarter
#' (`tail.flag = 0.25`) of the taxon range?
#' @param x A stratigraphic occurrence dataset
#' @param rank The column name in x containing the taxon
#' names for which trimmed ranges will be calculated
#' @param srt A column name in x denoting the occurrence FADs
#' @param end A column name in x denoting the occurrence LADs
#' @param step A positive integer specifying the time window size
#' (i.e. the duration represented by each row in the output matrix)
#' @param density A positive numeric specifying the step size for
#' densifying records. This should ideally be smaller than step
#' @param top The percentage by which the top of the range will
#' be trimmed
#' @param bottom The percentage by which the bottom of the range
#' will be trimmed
#' @param tail.flag a numeric vector of proportions in the range
#' 0 > x > 1 which will be used to test for long tails
#' @param method The method for quantifying occurrence density. By
#' default both histogram and kernel density will be used
#' @return If the user specifies a specific method (e.g. method = "kernel"),
#' the returned value will be a data.frame containing the taxa as row names,
#' the original taxon ranges (FAD, LAD), their ranges as trimmed by the
#' specified value (default FAD95, LAD95), and the tail status (0 = none, 1 = tail)
#' at the user-specified tail proportions. If method is not specified, the result
#' will be a list of 2 data.frames, one for each method
#' @import pbapply
#' @import data.table
#' @export
#' @source Pacman procedure modified from https://rdrr.io/github/plannapus/CONOP9companion/src/R/pacman.R.
#' @references Lazarus et al (2012) Paleobiology
#' @examples
#' # load dataset
#' data("brachios")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # run pacmacro
#' pacm <- pacmacro_ranges(brachios, tail.flag = c(0.3, 0.35, 0.4),
#'                         rank = "genus", srt = "max_ma", end = "min_ma")

pacmacro_ranges <- function (x, rank = "genus", srt = "max_ma", end = "min_ma",
                           step = 1, density = 0.1, top = 5, bottom = 5, tail.flag = 0.35, method = c("histogram", "kernel")) {

  if (!is.data.frame(x)) {
    stop("x must be a dataframe")
  }
  if (!all(c(rank, srt, end) %in% colnames(x))) {
    stop("One or more of rank, srt or end are not colnames in data")
  }
  if (!is.numeric(step) | !is.numeric(density)) {
    stop("Step and density must be numeric")
  }
  if (length(step) > 1) {
    stop("Step must be a single positive integer")
  }
  if (step < 1) {
    stop("Step must be a single positive integer")
  }
  if (length(density) > 1) {
    stop("Density must be a single positive integer")
  }
  if (density >= step) {
    warning("Density should ideally be smaller than step")
  }
  if (!is.numeric(x[, srt]) | !is.numeric(x[, end])) {
    stop("Columns srt and end must be numeric")
  }
  if (any(x[, srt] < x[, end])) {
    stop("One or more maximum ages are smaller than their corresponding minimum ages")
  }
  if (!all(method %in% c("histogram", "kernel"))) {
    stop("Method must be one or both of histogram or kernel")
  }
  . <- NULL
  x <- x[, c(rank, srt, end)]
  x <- x[!is.na(x[, rank]), ]
  #x[, srt] <- x[, srt] + density
  start_v <- ceiling(max(x[[srt]]) + density)
  end_v <- floor(min(x[[end]]))
  x <- data.table::as.data.table(x)
  data.table::setkeyv(x, rank)

  to_do <- stats::na.omit(unique(x[[rank]]))
  test <- pbsapply(to_do, simplify = FALSE, function(y) {

    upr <- x[.(y), c("max_ma", "min_ma")]
    occ_rng <- c(max(upr[[1]]), min(upr[[2]]))
    # small constant to make sure density method works
    upr[[1]] <- upr[[1]] + density
    upr <- as.vector(unlist(apply(upr, 1, function(z) {
      seq(from = z[1], to = z[2], by = -density)
    })))
    start_range <- ceiling((max(upr)))
    end_range <- floor((min(upr)))
    bins <- seq(from = end_range, to = start_range, by = step)
    pre_base <- NULL
    if (start_range != start_v) {
      pre_base <- rep(0, length(seq(from = start_v, to = start_range + step, by = -step)))
    }
    post_top <- NULL
    if (end_range != end_v) {
      post_top <- rep(0, length(seq(from = end_range - step, to = end_v, by = -step)))
    }

    # histogram
    bounds_h <- NULL
    if ("histogram" %in% method | "kernel" %in% method) {

      # pacman count
      upr_h0 <- rev(hist(upr, breaks = bins, plot = FALSE)$counts)
      upr_h <- c(pre_base, upr_h0, post_top)
      names(upr_h) <- seq(from = start_v, to = end_v + step, by = -step)
      nb_top <- sum(upr_h) * top / 100
      nb_bottom <- sum(upr_h) * bottom / 100
      #outlier_h <- rep(0, length(upr_h))
      trimmed_h <- upr_h
      trimmed_h[!(cumsum(trimmed_h) >= nb_top & rev(cumsum(rev(trimmed_h))) >= nb_bottom)] <- 0
      #outlier_h <- outlier_h - trimmed_h

      # skew count
      #h1 <- which(diff(cumsum(upr_h)) != 0)
      #h1 <- c(max(upr), min(upr))
      h1 <- occ_rng
      h2 <- which(diff(cumsum(trimmed_h)) != 0)
      h2 <- as.numeric(names(c(h2[1], h2[length(h2)]) + 1))
      if(h2[1] > h1[1] | is.infinite(h2[1])) {
        h2[1] <- h1[1]
      }
      if(h2[2] < h1[2] | is.infinite(h2[1])) {
        h2[2] <- h1[2]
      }
      if(h1[1] - h2[1] <= step) {h2[1] <- h1[1]}
      if(h2[2] - h1[2] <= step) {h2[2] <- h1[2]}
      bounds_h <- c(h1, h2)
    }

    bounds_k <- NULL
    if ("kernel" %in% method) {
      if (length(bins) < 3) {
        bounds_h <- bounds_k
      } else {
        den <- density(upr, n = (length(bins) - 1),
                       from = end_range, to = start_range)
        if (sum(den$y) == 0) {
          den$y <- rep(1/length(bins - 1), times = length(bins) - 1)
        }
        den$y[which(upr_h0 == 0)] <- 0
        # skew count
        #k1 <- c(max(upr), min(upr))
        k1 <- occ_rng
        k2 <- as.vector(c(ceiling(quantile_coef_density_BMS(den,
                                                             probs = 1 - top/100)), floor(quantile_coef_density_BMS(den,
                                                                                                                     probs = bottom/100))))
        if(k2[1] > k1[1] | is.infinite(k2[1])) {
          k2[1] <- k1[1]
        }
        if(k2[2] < k1[2] | is.infinite(k2[1])) {
          k2[2] <- k1[2]
        }
        if(k1[1] - k2[1] <= step) {k2[1] <- k1[1]}
        if(k2[2] - k1[2] <= step) {k2[2] <- k1[2]}
        bounds_k <- c(k1, k2)
        # pacman count
        #upr_k <- den$y
        #outlier <- upr_k
        #trimmed <- upr_k
        #trimmed[!(den$x < k2[2] & den$x > k2[1])] <- 0
        #trimmed <- c(pre_base, rev(trimmed), post_top)
        #outlier[!(den$x > k2[2] & den$x > k2[1])] <- 0
        #outlier <- c(pre_base, rev(outlier), post_top)
      }
    }
    bounds <- c(bounds_h, bounds_k)
  })

  # output
  test2 <- do.call(rbind, test)
  if (length(method) == 1) {
    if(method == "histogram") {
      out <- test2[,1:4]
    }
    if(method == "kernel") {
      out <- test2[,5:8]
    }
    tails <- do.call(cbind, lapply(tail.flag, function(z) {
      dur <- abs(out[,1] - out[,2])
      fdist <- out[,1] - out[,3]
      ldist <- out[,4] - out[,2]
      tdist <- fdist + ldist
      tail_logical <- as.numeric((tdist / dur) >= z)
    }))
    out <- cbind.data.frame(out, tails)
    tnam <- paste0("tflag", tail.flag)
    colnames(out) <- c("FAD", "LAD", paste0("FAD", 100 - bottom), paste0("LAD", 100 - top), tnam)
  }
  if (length(method) == 2) {
    o1 <- test2[,1:4]
    tails <- do.call(cbind, lapply(tail.flag, function(z) {
      dur <- abs(o1[,1] - o1[,2])
      fdist <- o1[,1] - o1[,3]
      ldist <- o1[,4] - o1[,2]
      tdist <- fdist + ldist
      tail_logical <- as.numeric((tdist / dur) >= z)
    }))
    o1 <- cbind.data.frame(o1, tails)
    o2 <- test2[,5:8]
    tails <- do.call(cbind, lapply(tail.flag, function(z) {
      dur <- abs(o2[,1] - o2[,2])
      fdist <- o2[,1] - o2[,3]
      ldist <- o2[,4] - o2[,2]
      tdist <- fdist + ldist
      tail_logical <- as.numeric((tdist / dur) >= z)
    }))
    o2 <- cbind.data.frame(o2, tails)
    tnam <- paste0("tflag", tail.flag)
    colnames(o1) <- colnames(o2) <- c("FAD", "LAD", paste0("FAD", 100 - bottom), paste0("LAD", 100 - top), tnam)
    out <- list(o1, o2)
    names(out) <- c("histogram", "kdensity")
  }
  return(out)
}
