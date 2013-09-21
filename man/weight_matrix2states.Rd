\name{weight_matrix2states}
\alias{weight_matrix2states}
\title{Returns unique state assignment from a (row-wise) weight matrix}
\usage{
weight_matrix2states(weight_matrix, rule = "argmax")
}
\arguments{
  \item{weight_matrix}{an \eqn{N \times K} matrix}

  \item{rule}{how do we choose the state given the weight
  matrix. \code{c("argmax", "sample")}.}
}
\description{
  Converts a probabilistic cluster assignment to a unique
  cluster assignment using the \describe{
  \item{\code{'argmax'} rule:}{state of row \eqn{i} is
  assigned as the position of the maximum in that row (ties
  are broken at random).} \item{\code{'sample'} rule}{state
  of row \eqn{i} is sampled from the discrete distribution
  where probabilities equal the weight vector in row
  \eqn{i}} }
}
\examples{
WW <- matrix(runif(12), ncol = 3)
WW <- normalize(WW)
WW
weight_matrix2states(WW)
weight_matrix2states(WW, "sample")
# another 'sample' is in general different from previous conversion unless WW
# is a 0/1 matrix
weight_matrix2states(WW, "sample")
}
\seealso{
  \code{\link{states2weight_matrix}}
}
\keyword{array}
\keyword{manip}
