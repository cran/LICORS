\name{LC-utils}
\alias{LC-utils}
\alias{LC_coords2control_settings}
\alias{plot.LC}
\alias{summary.LC}
\title{Utilities for LC class}
\usage{
\method{plot}{LC}(x, cex.axis = 2, cex.lab = 2, ...)

\method{summary}{LC}(object, verbose = TRUE, ...)

LC_coords2control_settings(LC_coords)
}
\arguments{
  \item{x}{an object of class \code{"LC"} (see
  \code{\link{setup_LC_geometry}})}

  \item{cex.axis}{The magnification to be used for axis
  annotation relative to the current setting of
  \code{cex}.}

  \item{cex.lab}{The magnification to be used for x and y
  labels relative to the current setting of \code{cex}.}

  \item{...}{optional arguments passed to \code{plot}.}

  \item{object}{an object of class \code{"LC"}}

  \item{verbose}{logical; if \code{TRUE} LC information is
  printed in the console}

  \item{LC_coords}{template of a light cone (with respect
  to origin)}
}
\description{
  The \code{"LC"} class is the core property of the LICORS
  model as it specifies the spatio-temporal neighborhood of
  the past and future light cone. The function
  \code{\link{setup_LC_geometry}} generates an instance of
  the \code{"LC"} class.

  \code{plot.LC} plots LCs of \eqn{(1+1)D} and \eqn{(2+1)D}
  systems with respect to the origin \eqn{(\mathbf{r}, t) =
  (\boldsymbol 0, 0)}. This is especially useful for a
  quick check if the LC geometry specified by
  \code{\link{setup_LC_geometry}} is really the intended
  one.

  \code{summary.LC} prints a summary of the LC geometry.
  Returns (invisible) the summary matrix.

  \code{LC_coords2control_setting} computes auxiliary
  measures given the LC geometry such as horizon,
  spatial/temporal extension, etc. This function should not
  be called by the user directly; only by
  \code{\link{get_spacetime_iterator}}.
}
\examples{
aa <- setup_LC_geometry(horizon = list(PLC = 2, FLC = 1), speed = 1, space_dim = 1, 
    shape = "cone")
plot(aa)
bb <- setup_LC_geometry(horizon = list(PLC = 2, FLC = 1), speed = 1, space_dim = 1, 
    shape = "revcone")
plot(bb)
aa <- setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed = 1, space_dim = 1, 
    shape = "cone")
summary(aa)
aa <- setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed = 1, space_dim = 1, 
    shape = "cone")
LC_coords2control_settings(aa$coords)
}
\seealso{
  \code{\link{compute_LC_coords}}
}
\keyword{hplot}
\keyword{method}
\keyword{models}
\keyword{print}
