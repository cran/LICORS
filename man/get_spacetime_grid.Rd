\name{get_spacetime_grid}
\alias{get_spacetime_grid}
\title{Get an iterator over the space-time coordinates of the field.}
\usage{
get_spacetime_grid(dim, LC.coordinates)
}
\arguments{
  \item{dim}{dimension of the original field (first
  dimension is time; rest is space)}

  \item{LC.coordinates}{template of the LC coordinates}
}
\description{
  This function returns a matrix of space-time coordinates
  of the field. Both for the whole field as well as the
  truncated field (without the margin)
}
\examples{
AA <- matrix(rnorm(200), ncol = 10)
LC.geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 3, FLC = 0), shape = "cone")
bb <- get_spacetime_grid(dim(AA), LC.geom$coordinates)
}
\seealso{
  \code{\link{compute_LC_coordinates}},
  \code{\link{setup_LC_geometry}}
}
