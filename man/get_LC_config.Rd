\name{get_LC_config}
\alias{get_LC_config}
\title{Get configuration of a light cone (LC)}
\usage{
get_LC_config(coord, field, LC.coordinates)
}
\arguments{
  \item{coord}{space-time coordinate \eqn{(\mathbf{r}, t)}}

  \item{field}{spatio-temporal field; either a matrix or a
  3-dimensional array with time \eqn{t} as the first coord,
  and the spatial coords in order. Make sure to see also
  \code{\link{compute_LC_coordinates}} for correct
  formatting.}

  \item{LC.coordinates}{template coords for the LC}
}
\description{
  \code{get_LC_config} obtains the PLC or FLC at a
  particular \eqn{(\mathbf{r}, t)} from a \eqn{(N+1)D}
  field based on the LC template from
  \code{\link{compute_LC_coordinates}} (or
  \code{\link{setup_LC_geometry}}).
}
\examples{
AA <- matrix(rnorm(40), ncol = 5)
image2(AA)
LCind <- compute_LC_coordinates(speed = 1, horizon = 1, shape = "cone")
AA
get_LC_config(cbind(5, 2), AA, LCind)
# a time series example
data(nhtemp)
xx <- c(nhtemp)
LCind <- compute_LC_coordinates(speed = 1, horizon = 4, shape = "cone", space.dim = 0)
cc <- get_LC_config(6, xx, LCind)
}
\seealso{
  \code{\link{compute_LC_coordinates}}
}
\keyword{method}
\keyword{utilities}
