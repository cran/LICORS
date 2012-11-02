\name{compute_margin_coords}
\alias{compute_margin_coords}
\title{Get LC configuration from a (N+1)D field}
\usage{
compute_margin_coords(dim = NULL, LC_coords = NULL)
}
\arguments{
  \item{dim}{a vector with the dimensions of the field
  (time, space1, space2, ..., spaceN)}

  \item{LC_coords}{template of the LC coords}
}
\description{
  \code{compute_margin_coords} computes the coordinates
  (boundary) of the margin of the field.
}
\examples{
LC_geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 3, FLC = 0), shape = "cone")

data(contCA00)

aa <- compute_margin_coords(dim(contCA00$observed), LC_geom$coords)
aa
}
\seealso{
  \code{\link{compute_LC_coords}}
}
\keyword{manip}
