\name{compute_LC_coords}
\alias{compute_LC_coords}
\title{Computes coordinates of PLC and FLC relative to origin}
\usage{
compute_LC_coords(horizon = 1, speed = 1, space_dim = 1, type = c("PLC", "FLC"), 
    shape = c("cone", "tube", "revcone"))
}
\arguments{
  \item{horizon}{integer; horizon for the PLC or FLC}

  \item{speed}{speed of propagation}

  \item{space_dim}{maximum value}

  \item{type}{\code{"PLC"} or \code{"FLC"}}

  \item{shape}{shape of light cone: \code{'cone'},
  \code{'tube'}, or \code{'revcone'}.}
}
\description{
  Computes the space-time coordinates of PLC and FLC given
  control settings relative to the origin \eqn{(\mathbf{r},
  t) = (\boldsymbol 0, 0)}.

  Since these coordinates do not change for different
  space-time positions, they can be computed once before
  getting the LC configurations for the entire field and
  then used in each call by array maskexing in
  \code{\link{get_LC_config}}.
}
\examples{
plot(compute_LC_coords(speed = 1, horizon = 4), xlim = c(-4, 2), pch = "-", cex = 2, 
    col = 2, xlab = "Time", ylab = "Space")
points(compute_LC_coords(speed = 1, horizon = 2, type = "FLC"), pch = "+", cex = 2, 
    col = "blue")

plot(compute_LC_coords(speed = 1, horizon = 4, shape = "tube", type = "FLC"))
plot(compute_LC_coords(speed = 1, horizon = 4, shape = "revcone", type = "PLC"))
}
\seealso{
  \code{\link{get_LC_config}}
  \code{\link{setup_LC_geometry}} \code{\link{summary.LC}}
  \code{\link{plot.LC}}
}
\keyword{manip}
