\name{states2weight_matrix}
\alias{states2weight_matrix}
\title{Converts label vector to 0/1 mixture weight matrix}
\usage{
states2weight_matrix(state_vector = NULL, nstates_total = NULL)
}
\arguments{
  \item{state_vector}{a vector of length \eqn{N} with the
  state labels}

  \item{nstates_total}{total number of states. If
  \code{NULL}, then the maximum of \code{state_vector} is
  chosen}
}
\description{
  Converts unique cluster assignment stored in a length
  \eqn{N} label vector into a \eqn{N \times K} Boolean
  matrix of mixture weights.
}
\examples{
ss <- sample.int(5, 10, replace = TRUE)
WW <- states2weight_matrix(ss)

image2(WW, col = "RdBu", xlab = "States", ylab = "Samples", axes = FALSE)
}
\seealso{
  \code{\link{weight_matrix2states}}
}
\keyword{array}
\keyword{manip}
