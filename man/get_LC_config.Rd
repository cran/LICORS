\name{get_LC_config}
\alias{get_LC_config}
\title{Get configuration of a light cone (LC)}
\usage{
get_LC_config(coord, field = NULL, LC_coords = NULL)
}
\arguments{
  \item{coord}{space-time coordinate \eqn{(\mathbf{r}, t)}}

  \item{field}{spatio-temporal field; either a matrix or a
  3-dimensional array with time \eqn{t} as the first coord,
  and the spatial coords in order. Make sure to see also
  \code{\link{compute_LC_coords}} for correct formatting.}

  \item{LC_coords}{template coords for the LC}
}
\description{
  \code{get_LC_config} obtains the PLC or FLC at a
  particular \eqn{(\mathbf{r}, t)} from a \eqn{(N+1)D}
  field based on the LC template from
  \code{\link{compute_LC_coords}} (or
  \code{\link{setup_LC_geometry}}).
}
\examples{
AA <- matrix(rnorm(40), ncol = 5)
image2(AA)
LCind <- compute_LC_coords(speed = 1, horizon = 1, shape = "cone")
AA
get_LC_config(c(5, 2), AA, LCind)
}
\seealso{
  \code{\link{compute_LC_coords}}
}
\keyword{method}
\keyword{utilities}
