% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_coef_density_BMS.R
\name{quantile_coef_density_BMS}
\alias{quantile_coef_density_BMS}
\title{quantile_coef_density_BMS}
\source{
static rip from BMS package
}
\usage{
quantile_coef_density_BMS(
  x,
  probs = seq(0.25, 0.75, 0.25),
  names = TRUE,
  normalize = TRUE,
  ...
)
}
\arguments{
\item{x}{a object of class pred.density, coef.density, density, or a list of densities}

\item{probs}{numeric vector of probabilities with values in range 0 - 1. Elements very close to the boundaries return Inf or -Inf}

\item{names}{logical; if TRUE, the result has a names attribute, resp. a rownames and colnames attributes. Set to FALSE for speedup with many probs}

\item{normalize}{logical if TRUE then the values in x$y are multiplied with a factor such that their integral is equal to one}

\item{...}{further arguments passed to or from other methods.}
}
\value{
If x is of class density (or a list with exactly one element), a vector with quantiles.
If x is a list of densities with more than one element (e.g. as resulting from pred.density or coef.density),
then the output is a matrix of quantiles, with each matrix row corresponding to the respective density.
}
\description{
Static rip of the quantile.coef.density
function and relevant internals from the
BMS package as the package is archived.
}
