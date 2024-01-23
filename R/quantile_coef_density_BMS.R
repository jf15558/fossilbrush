#' quantile_coef_density_BMS
#'
#' Static rip of the quantile.coef.density
#' function and relevant internals from the
#' BMS package as the package is archived.
#' @param x a object of class pred.density, coef.density, density, or a list of densities
#' @param probs numeric vector of probabilities with values in range 0 - 1. Elements very close to the boundaries return Inf or -Inf
#' @param names logical; if TRUE, the result has a names attribute, resp. a rownames and colnames attributes. Set to FALSE for speedup with many probs
#' @param normalize logical if TRUE then the values in x$y are multiplied with a factor such that their integral is equal to one
#' @param ... further arguments passed to or from other methods.
#' @return If x is of class density (or a list with exactly one element), a vector with quantiles.
#' If x is a list of densities with more than one element (e.g. as resulting from pred.density or coef.density),
#' then the output is a matrix of quantiles, with each matrix row corresponding to the respective density.
#' @source static rip from BMS package

quantile_coef_density_BMS <- function (x, probs = seq(0.25, 0.75, 0.25), names = TRUE, normalize = TRUE, ...) {

  # quick inline function
  my.quantile.density = function(x, probs, names, normalize = normalize, ...) {
    ycs = (cumsum(x$y) - (x$y - x$y[[1]])/2) * diff(x$x[1:2])
    if(normalize){
      ycs = ycs/(ycs[[length(ycs)]])
    }
    xin = x$x
    maxi = length(ycs)
    qqs = sapply(as.list(probs), function(qu) {
      iii = sum(ycs <= qu)
      if(iii == maxi) {
        return(Inf)
      } else if(iii == 0L) {
        return(-Inf)
      } else {
        return(xin[[iii + 1]] + ((ycs[[iii + 1]] - qu)/(ycs[[iii + 1]] - ycs[[iii]])) * (xin[[iii]] - xin[[iii + 1]]))
      }
    })
    if(as.logical(names)) {
      names(qqs) = paste(format(100 * probs, trim = TRUE, digits = max(2L, getOption("digits"))), "%", sep = "")
    }
    return(qqs)
  }

  probs = as.vector(probs)
  if(is.element("density", class(x))) {
    return(my.quantile.density(x = x, probs = probs, names = names, normalize = normalize))
  }
  if(!all(sapply(x, function(dd) is.element("density", class(dd))))) {
    stop("x needs to be a density or list of densities")
  }
  if(length(x) == 1L) {
    return(my.quantile.density(x = x[[1]], probs = probs, names = names, normalize = normalize))
  }

  qout = sapply(x, my.quantile.density, probs = probs, names = FALSE, normalize = normalize)
  if(!is.matrix(qout)) {
    if (length(probs) > 1) {
      return(qout)
    }
    qout = as.matrix(qout)
  } else {
    qout = t(qout)
  }
  if(as.logical(names)) {
    colnames(qout) = paste(format(100 * probs, trim = TRUE, digits = max(2L, getOption("digits"))), "%", sep = "")
  }
  if(is.matrix(quout) && as.logical(names)) {
    rownames(quout) <- sapply(x, function(lx) lx[["data.name"]])
  }
  return(quout)
}
