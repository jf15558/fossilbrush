#' densify
#'
#' Function to create a matrix of occurrence record densities
#' through geological time from an occurrence dataset. Each column
#' represents a taxon. Each row represents a user defined window of
#' time, with the first row starting at the oldest FAD in the
#' dataset and spanning to the youngest LAD stepwise by the user
#' defined window (default of 1 Ma). Occurrence records are
#' densified by generating a vector of time points from occurrence
#' FAD to occurrence LAD (default step of 0.1 Ma), then tallied in
#' two ways. The first way is a simple histogram count of
#' points-per-window, with the same number of histogram bins as time
#' steps between the overall taxon FAD and LAD. The second way is a
#' kernel density estimate, using a Gaussian kernel with a equally
#' spaced estimatopms equal to the number of timesteps between the
#' overall taxon FAD and LAD
#' @param x An occurrence dataset
#' @param rank The column name in x containing the taxon names
#' for which densified columns will be generated
#' @param srt A column name in x denoting the occurrence FADs
#' @param end A column name in x denoting the occurrence LADs
#' @param step A positive integer specifying the time window size
#' (i.e. the duration represented by each row in the output matrix)
#' @param density A positive numeric specifying the step size for
#' densifying records. This should ideally be smaller than step
#' @param method The method for quantifying occurrence density. By
#' default both histogram and kernel density will be used
#' @param ... additional arguments passed to @seealso density
#' @param verbose A logical determining if function progress
#' should be reported
#' @return A list of two sparse matrices, the first containing
#' the histogram counts, the second the kernel density estimates
#' @import pbapply
#' @import data.table
#' @importFrom stats na.omit
#' @importFrom methods as slot
#' @importClassesFrom Matrix sparseVector
#' @export

densify <- function(x, rank = "genus", srt = "max_ma", end = "min_ma", step = 1, density = 0.1,
                    method = c("histogram", "kernel"), ..., verbose = TRUE) {

  if(!is.data.frame(x)) {
    stop("x must be a dataframe")
  }
  if(!all(c(rank, srt, end) %in% colnames(x))) {
    stop("One or more of rank, srt or end are not colnames in data")
  }
  if(class(step) != "numeric" | class(density) != "numeric") {
    stop("Step and density must be numeric")
  }
  if(length(step) > 1) {
    stop("Step must be a single positive integer")
  }
  if(step < 1) {
    stop("Step must be a single positive integer")
  }
  if(length(density) > 1) {
    stop("Density must be a single positive integer")
  }
  if(density >= step) {
    warning("Density should ideally be smaller than step")
  }
  if(class(x[,srt]) != "numeric" | class(x[,end]) != "numeric") {
    stop("Columns srt and end must be numeric")
  }
  if(any(x[,srt] < x[,end])) {
    stop("One or more maximum ages are smaller than their corresponding minimum ages")
  }
  if(!all(method %in% c("histogram", "kernel"))) {
    stop("Method must be one or both of histogram or kernel")
  }
  if(!is.logical(verbose) | length(verbose) != 1) {
    stop("verbose should be a logical of length 1")
  }
  if(!verbose) {
    baseopt <- getOption("pboptions")
    opb <- pboptions(type = "none")
  }

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

  # cut down to columns
  x <- x[,c(rank, srt, end)]
  x <- x[!is.na(x[,rank]),]
  # add a very small constant to deal with zero range ages
  x[,srt] <- x[,srt] + density
  # set the upper and lower range for the vectors
  start_v <- ceiling(max(x[[srt]]))
  end_v <- floor(min(x[[end]]))
  x <- data.table::as.data.table(x)
  data.table::setkeyv(x, rank)

  # get the unique names for which to get densified ranges
  to_do <- stats::na.omit(unique(x[[rank]]))

  # create the densified occurrence record
  test <- pbsapply(to_do, simplify = FALSE, function(y) {

    # get all occurrences of the taxon (uses data.table)
    upr <- x[.(y), c("max_ma", "min_ma")]
    # sequence from each FAD-LAD pair by density
    upr <- as.vector(unlist(apply(upr, 1, function(z) {seq(from = z[1], to = z[2], by = -density)})))
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
    upr_k <- NULL
    if("kernel" %in% method) {

      # if only 1 bin, put all density in that bin
      if(length(bins) < 3) {
        upr_k <- c(pre_base, 1, post_top)
      } else {
        # otherwise do density across bins
        den <- rev(density(upr, n = (length(bins) - 1), from = end_range, to = start_range, ...)$y)
        # if statement to catch very rare 0 density estimates (unclear why, but only 6 cases out of 66000)
        if(sum(den) == 0) {
          den <- rep(1 / length(bins - 1), times = length(bins) - 1)
        }
        upr_k <- c(pre_base, den, post_top)
        upr_k[which(upr_h == 0)] <- 0
      }
    }
    # bind and store as sparse vector
    upr2 <- c(upr_h, upr_k)
    upr2 <- as(upr2, "sparseVector")
  })
  test2 <- do.call(sv_cbind, test)

  # format output
  if(length(method) == 1) {
    o1 <- test2[((nrow(test2) / 2) + 1):nrow(test2),]
    dimnames(o1) <- list(seq(from = start_v, to = end_v + step, by = -step), names(test))
    out <- list(o1)
    if(method == "kernel") {method <- "kdensity"}
    names(out) <- method
  }
  if(length(method) == 2) {
    o1 <- test2[1:(nrow(test2) / 2),]
    o2 <- test2[((nrow(test2) / 2) + 1):nrow(test2),]
    dimnames(o2) <- dimnames(o1) <- list(seq(from = start_v, to = end_v + step, by = -step), names(test))
    out <- list(o1, o2)
    names(out) <- c("histogram", "kdensity")
  }
  if(!verbose) {opt <- pboptions(baseopt)}
  return(out)
}
