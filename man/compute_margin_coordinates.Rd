\name{compute_margin_coordinates}
\alias{compute_margin_coordinates}
\title{Get LC configuration from a (N+1)D field}
\usage{
compute_margin_coordinates(dim, LC.coordinates)
}
\arguments{
  \item{dim}{a vector with the dimensions of the field
  (time, space1, space2, ..., spaceN)}

  \item{LC.coordinates}{template of the LC coordinates}
}
\description{
  \code{compute_margin_coordinates} computes the
  coordinates (boundary) of the margin of the field.
}
\examples{
LC_geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 3, FLC = 0), shape = "cone")

data(contCA00)

aa <- compute_margin_coordinates(dim(contCA00$observed), LC_geom$coordinates)
aa
}
\seealso{
  \code{\link{compute_LC_coordinates}}
}
\keyword{manip}
