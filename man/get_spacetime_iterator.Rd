\name{get_spacetime_iterator}
\alias{get_spacetime_iterator}
\title{Get an iterator over the space-time coordinates of the field}
\usage{
get_spacetime_iterator(dim = NULL, LC_coords = NULL)
}
\arguments{
  \item{dim}{dimension of the original field (first
  dimension is time; rest is space)}

  \item{LC_coords}{template of the LC coords}
}
\description{
  This function returns iterators over space-time
  coordinates of the field. Both for the whole field as
  well as the truncated field (without the margin)
}
\examples{
AA <- matrix(rnorm(200), ncol = 10)
LC_geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 3, FLC = 0), shape = "cone")
bb <- get_spacetime_iterator(dim(AA), LC_geom$coords)
}
\seealso{
  \code{\link{compute_LC_coords}},
  \code{\link{setup_LC_geometry}}
}
\keyword{iteration}
