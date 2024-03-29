% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chrono_scale.R
\name{chrono_scale}
\alias{chrono_scale}
\title{chrono_scale}
\usage{
chrono_scale(
  x,
  tscale = "GTS_2020",
  srt = "early_interval",
  end = "late_interval",
  max_ma = NULL,
  min_ma = NULL,
  verbose = TRUE
)
}
\arguments{
\item{x}{A data.frame containing, minimally two columns
corresponding respectively to the first and last intervals
of the data. Values should only be present in the second column
where the minimum age interval for a row is different to the
maximum age interval. Otherwise the values should be NA and the
ages returned will be based on the interval specified in the
first column, in line with PBDB formatting.}

\item{tscale}{A character string specifying one of the inbuilt
chronostratigraphic timescales (currently GTS 2020 only) or a
data.frame supplied by the user. If the latter, this must
contain columns named 'Interval', 'FAD', 'LAD', specifying the
interval names to be matched and their lower and upper age in Ma}

\item{srt}{A character of length 1 specifing the column name
of the first interval field in x}

\item{end}{A character of length 1 specifing the column name
of the last interval field in x}

\item{max_ma}{If not NULL, a character of length 1 specifing
the column name of the original numeric maximum age field in
x, to be used as fall back values if interval names cannot
all be matched}

\item{min_ma}{If not NULL, a character of length 1 specifing
the column name of the original numeric minimum age field in
x, to be used as fall back values if interval names cannot
all be matched}

\item{verbose}{A logical indicating if warning messages should be
displayed or otherwise}
}
\value{
The dataframe, x, with two additional columns containing
the revised first and last numeric ages of the data, with column
names GTS_FAD and GTS_LAD respectively
}
\description{
Convenience function to apply user-specified chronostratigraphy
to fossil datasets. The function relies on a lookup table
generated based on the named intervals in the PBDB in early
2021. First and last interval names in the supplied dataset
are matched against this lookup table, by default using
'get("GTS_2020)", to get GTS_2020 numeric ages. If the
dataset contains intervals which are not present in the
lookup table, they will not be matched and the user will
be warned. To get around this possibility, the user can
also supply the original numeric ages which will be used
as default ages if an interval cannot be matched, to ensure
that the returned vectors of numeric ages do not contain NAs.
}
\examples{
# example dataset
data("brachios")
# add GTS_2020 dates
brachios <- chrono_scale(brachios, srt = "early_interval", end = "late_interval",
                          max_ma = "max_ma", min_ma = "min_ma")
}
