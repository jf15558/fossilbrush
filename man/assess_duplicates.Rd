% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assess_duplicates.R
\name{assess_duplicates}
\alias{assess_duplicates}
\title{assess_duplicates}
\usage{
assess_duplicates(
  x,
  node,
  mode = c("frequency", "completeness"),
  jump = 3,
  plot = FALSE
)
}
\arguments{
\item{x}{A tgraph object}

\item{node}{character vector of elements with multiple higher
classifications in x, or a tvertseq object with those same
elements as focal}

\item{mode}{The rule to be used in selecting between multiple
higher classifications. It is possible for the most complete
pathway to also be the most frequent}

\item{jump}{The maximum number of levels between the point of
divergence and the point of reunion (if present) for a given
path, below which the divergence will be taken as conflicting}

\item{plot}{A logical speciying if the divergent paths should
be plotted}
}
\value{
A list with as many items as elements with multiple
classifications, each recording the assessment for a given
element
}
\description{
Function to assess and resolve elements with multiple higher
classifications in a tgraph object. Assessment is performed
based on the topology of the graph they form. Linear paths
(i.e. two totally separate paths diverging from the a shared
node), rings (divergent paths which only reunite at the
highest rank in the tgraph) or more than two divergent paths
are treated as distinct. If not any of these cases, the
distance between the focal element and the reunion of the
divergent paths, along with their subtopologies are assessed
and a consensus or preferred path based on the frequency of
each path in the tgraph or their completeness returned, or
the element judged as having multiple distinct
classifications
}
