\name{sparsify_weights}
\alias{sparsify_weights}
\title{Sparsify weights}
\usage{
sparsify_weights(weight.matrix.proposed, weight.matrix.current = NULL, 
    penalty = "entropy", lambda = 0)
}
\arguments{
  \item{weight.matrix.proposed}{\eqn{N \times K} weight
  matrix}

  \item{weight.matrix.current}{\eqn{N \times K} weight
  matrix}

  \item{penalty}{type of penalty: \code{c("entropy",
  "1-Lq", "lognorm")}. Default: \code{"entropy"}}

  \item{lambda}{penalization parameter: larger
  \code{lambda} gives sparser mixture weights}
}
\description{
  This function makes weights of a mixture model more
  sparse using gradient based penalty methods.
}
\examples{
WW <- matrix(c(rexp(10, 1/10), runif(10)), ncol = 5, byrow = FALSE)
WW <- normalize(WW)
WW_sparse <- sparsify_weights(WW, lambda = 0.1)
WW_more_sparse <- sparsify_weights(WW, lambda = 0.5)
compute_mixture_penalty(WW)
compute_mixture_penalty(WW_sparse)
compute_mixture_penalty(WW_more_sparse)
}
\seealso{
  \code{\link{compute_mixture_penalty}},
  \code{\link{mixed_LICORS}}
}
\keyword{array}
\keyword{manip}
