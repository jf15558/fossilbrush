#' add_itp
#'
#' Function to add detected peaks using the output of
#' @seealso threshold_ranges. This function should be used
#' to add information to an existing plot from @seealso densify,
#' ensuring that the same taxon name is being used
#' @param x The list output of @seealso threshold_ranges
#' @param taxon A character vector of length one, specifying
#' one of the taxon names in x to be plotted
#' @param legend.pos One of topleft, bottomleft, topright or
#' bottomright, or a vector of length two, giving the xy
#' coordinates of the legend. A convenience parameter so that
#' the plot detail can remain unobscured.
#' @param exit Restore base plotting parameters on function exit
#' (default as a requirement for CRAN). Can be set to false to allow
#' other elements to be aded to a plot
#' @return None, the detected peaks are added to an existing density plot
#' @export
#' @examples
#' # load dataset#
#' data("brachios")
#' # subsample brachios to make for a short example runtime
#' set.seed(1)
#' brachios <- brachios[sample(1:nrow(brachios), 1000),]
#' # densify ranges
#' dens <- densify(brachios)
#' # interpeak thresholding
#' itp <- threshold_ranges(brachios, win = 8, thresh = 10,
#'                         rank = "genus", srt = "max_ma", end = "min_ma")
#' # append the stratigraphically thresholded taxon names to the dataset
#  brachios$newgen <- itp$data
#' # plot the taxon, now identifying the peaks
#' plot_dprofile(dens, "Atrypa", exit = FALSE)
#' add_itp(itp, "Atrypa")

add_itp <- function(x, taxon, legend.pos = "topright", exit = TRUE) {

  if(exit) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  # check args
  if(!is.list(x)) {stop("x should be a list as outputted by 'threshold_peaks'")}
  if(length(x) > 4) {stop("x should be a list as outputted by 'threshold_peaks'")}
  if(!all(names(x) %in% c("data", "matrix", "peaks", "comparison"))) {stop("x should be a list as outputted by 'threshold_peaks'")}
  if(!is.vector(taxon) & length(taxon) == 1 & is.character(taxon)) {stop("taxon should be a character vector of length one")}
  if(!taxon %in% dimnames(x$matrix)[[2]]) {stop("taxon is not present in x")}
  if(is.numeric(legend.pos)) {
    legend.pos <- na.omit(legend.pos)
    if(length(legend.pos != 2)) {stop("If specifying legend.pos as a vector, this should be 2 positive numbers giving the xy position")}
    if(legend.pos[1] < 0 | legend.pos[2] < 0) {stop("If specifying legend.pos as a vector, this should be 2 positive numbers giving the xy position")}
  }
  if(is.character(legend.pos)) {
    legend.pos <- na.omit(legend.pos)
    if(length(legend.pos) != 1) {stop("If specifying legend.pos by character, this must be one of topleft, topright, bottomleft or bottomright")}
    if(!legend.pos %in% c("topleft", "topright", "bottomleft", "bottomright")) {stop("If specifying legend.pos by character, this must be one of topleft, topright, bottomleft or bottomright")}
  }

  # add lines
  for(i in seq(from = 1, to = length(x$peaks$profile_all[taxon][[1]]), by = 2)) {
    abline(v = names(x$peaks$profile_all[taxon][[1]])[i], col = "grey70")
    abline(v = names(x$peaks$profile_m[taxon][[1]])[i], lty = 2, col = "green")
    abline(v = names(x$peaks$profile_ms[taxon][[1]])[i], col = "green")
  }
  if(is.character(legend.pos)) {
    legend(x = legend.pos[1], y = legend.pos[2], lty = c(1, 2, 1), col = c("grey70", "green", "green"), legend = c("all", "sig. by mean", "sig. by mean + SD"), bty = "n")
  }
  if(is.numeric(legend.pos)) {
    legend(legend.pos, lty = c(1, 2, 1), col = c("grey70", "green", "green"), legend = c("all", "sig. by mean", "sig. by mean + SD"), bty = "n")
  }
}
