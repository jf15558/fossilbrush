#' plot_dprofile
#'
#' Function to plot density profiles of occurrences through
#' time using the output of @seealso densify.
#' @param x The list output of @seealso densify
#' @param taxon A character vector of length one, specifying
#' one of the taxon names in x to be plotted
#' @param exit Restore base plotting parameters on function exit
#' (default as a requirement for CRAN). Can be set to false to allow
#' other elements to be aded to a plot
#' @return NULL, the plotted density profile
#' @export
#' @examples
#' # load dataset
#' data("brachios")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # densify ranges
#' dens <- densify(brachios)
#  # plot an example taxon
#' plot_dprofile(dens, "Atrypa")

plot_dprofile <- function(x, taxon, exit = TRUE) {

  if(exit) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  if(!is.list(x)) {stop("x should be a list as outputted by 'densify'")}
  if(length(x) > 2) {stop("x should be a list as outputted by 'densify'")}
  if(!all(names(x) %in% c("kdensity", "histogram"))) {stop("x should be a list as outputted by 'densify'")}
  if(!is.vector(taxon) & length(taxon) == 1 & is.character(taxon)) {stop("taxon should be a character vector of length one")}
  if(!taxon %in% dimnames(x[[1]])[[2]]) {stop("taxon is not present in x")}

  lims <- x[[1]][,which(dimnames(x[[1]])[[2]] == taxon)]
  lims2 <- which(lims > 0)
  lims <- as.numeric(names(lims)[c(lims2[1], lims2[length(lims2)])])
  haxes <- TRUE
  hlab <- "Frequency"
  hlab2 <- ""
  nms <- names(x)
  par(xpd = FALSE)

  if(length(x) == 1) {
    par(mar = c(5, 5, 1, 1))
  } else {
    par(mar = c(5, 5, 1, 5))
  }
  if("kdensity" %in% names(x)) {
    plot(as.numeric(dimnames(x$kdensity)[[1]]),
         x$kdensity[,which(dimnames(x$kdensity)[[2]] == taxon)],
         type = "l", xlim = lims,
         ylab = "Density", xlab = "Time (Ma)")
    polygon(as.numeric(dimnames(x$kdensity)[[1]]),
            c(0, x$kdensity[2:(nrow(x$kdensity) - 1),which(dimnames(x$kdensity)[[2]] == taxon)], 0),
            col = "grey95")
  }
  if(length(x) == 2) {
    par(new = TRUE)
    haxes <- FALSE
    hlab <- ""
    hlab2 <- ""
  }
  if("histogram" %in% names(x)) {
    plot(as.numeric(dimnames(x$histogram)[[1]]),
         x$histogram[,which(dimnames(x$histogram)[[2]] == taxon)],
         type = "l", xlim = lims,
         col = "red", axes = haxes, xlab = hlab2, ylab = hlab, col.lab = "red")
  }
  if(length(x) == 2) {
    axis(4)
    mtext(side = 4, line = 3, "Frequency", col = "red")
  }
}
