\name{compute_mixture_penalty}
\alias{compute_mixture_penalty}
\title{Penalty of mixture weights}
\usage{
compute_mixture_penalty(weight_matrix, type = "entropy", q = 2, row_average = TRUE, 
    base = 2)
}
\arguments{
  \item{weight_matrix}{\eqn{N \times K} weight matrix}

  \item{type}{type of penalty: \code{c("entropy", "1-Lq",
  "lognorm")}. Default: \code{"entropy"}}

  \item{q}{exponent for \eqn{L_q} norm.}

  \item{row_average}{logical; if \code{TRUE} (default) then
  an average penalty over all rows will be returned (one
  single number); if \code{FALSE} a vector of length
  \eqn{N} will be returned.}

  \item{base}{logarithm base for the \code{"entropy"}
  penalty. Default: \code{base = 2}.  Any other real number
  is allowed; if \code{base = "nstates"} then it will
  internally assign it \code{base = ncol(weight_matrix)}.}
}
\description{
  Computes the penalty \eqn{\Omega(\mathbf{W})} of the
  weight matrix (or vector) for a mixture model.
}
\examples{
WW <- matrix(c(rexp(10, 1/10), runif(10)), ncol = 2, byrow = FALSE)
WW <- normalize(WW)
compute_mixture_penalty(WW, row_average = FALSE)
compute_mixture_penalty(WW, row_average = TRUE)  # default: average penalty
}
\seealso{
  \code{\link{compute_LICORS_loglik}}
  \code{\link{compute_NEC}}
}
\keyword{array}
\keyword{manip}
