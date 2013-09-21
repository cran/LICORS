\name{mixed_LICORS-utils}
\alias{mixed_LICORS-utils}
\alias{plot.mixed_LICORS}
\alias{summary.mixed_LICORS}
\title{Utilities for ``LICORS'' class}
\usage{
\method{plot}{mixed_LICORS}(x, type = "both", cex.axis = 1.5, cex.lab = 1.5, 
    cex.main = 2, line = 1.5, ...)

\method{summary}{mixed_LICORS}(object, ...)
}
\arguments{
  \item{x}{object of class \code{"mixed_LICORS"}}

  \item{type}{should only \code{"training"}, \code{"test"},
  or \code{"both"} be plotted. Default: \code{"both"}.}

  \item{cex.axis}{The magnification to be used for axis
  annotation relative to the current setting of
  \code{cex}.}

  \item{cex.lab}{The magnification to be used for x and y
  labels relative to the current setting of \code{cex}.}

  \item{cex.main}{The magnification to be used for main
  titles relative to the current setting of \code{cex}.}

  \item{line}{on which margin line should the labels be
  ploted, starting at 0 counting outwards (see also
  \code{\link[graphics]{mtext}}).}

  \item{...}{optional arguments passed to \code{plot} or to
  \code{summary}}

  \item{object}{object of class \code{"mixed_LICORS"}}
}
\description{
  The \code{"mixed_LICORS"} class is the output from the
  \code{\link{mixed_LICORS}} estimator.

  \code{plot.mixed_LICORS} gives a visual summary of the
  estimates such as marginal state probabilities,
  conditional state probabilities (= weight matrix),
  predictive state densities, trace plots for
  log-likelihood/loss/penalty.

  \code{summary.mixed_LICORS} prints out a summary of the
  estimated LICORS model.
}
\examples{
# see examples of LICORS-package see examples in LICORS-package
}
\keyword{hplot}
\keyword{model}
\keyword{nonparametric}
