\name{compute_NEC}
\alias{compute_NEC}
\title{Compute Negative Entropy Criterion (NEC)}
\usage{
compute_NEC(weight_matrix, loglik_1 = NULL, loglik_k = NULL)
}
\arguments{
  \item{weight_matrix}{\eqn{N \times K} weight matrix}

  \item{loglik_1}{baseline log-likelihood for \eqn{K=1}
  cluster model}

  \item{loglik_k}{log-likelihood for \eqn{K} cluster model}
}
\description{
  Computes the negative entropy criterion (NEC) to assess
  the number of clusters in a mixture model. See References
  for details.
}
\examples{
WW <- matrix(c(rexp(10, 1/10), runif(10)), ncol = 5, byrow = FALSE)
WW <- normalize(WW)
compute_NEC(WW, -2, -1)
}
\references{
  Christophe Biernacki, Gilles Celeux, and G\'erand
  Govaert(1999). ``An improvement of the NEC criterion for
  assessing the number of clusters in a mixture model''.
  Non-Linear Anal. 20, 3, 267-272.
}
\seealso{
  \code{\link{compute_mixture_penalty}}
}
\keyword{array}
\keyword{manip}
