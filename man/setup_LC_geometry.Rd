\name{setup_LC_geometry}
\alias{setup_LC_geometry}
\title{Setup light cone geometry}
\usage{
setup_LC_geometry(horizon = list(PLC = 1, FLC = 0), speed = 1, space_dim = 1, 
    shape = "cone")
}
\arguments{
  \item{horizon}{a list with \code{PLC} and \code{FLC}
  horizon}

  \item{speed}{speed of propagation}

  \item{space_dim}{dimension of the spatial grid. Eg.
  \code{2} if the data is a video ( = image sequences).}

  \item{shape}{shape of light cone: \code{'cone'},
  \code{'tube'}, or \code{'revcone'}.}
}
\value{
  A list of class \code{"LC"}.
}
\description{
  \code{setup_LC_geometry} sets up the light cone geometry
  for LICORS.
}
\examples{
aa <- setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed = 1, space_dim = 1, 
    shape = "cone")
aa
plot(aa)
summary(aa)
}
\seealso{
  \code{\link{LC-utils}}, \code{\link{compute_LC_coords}}
}
\keyword{manip}
