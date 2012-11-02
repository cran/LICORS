\name{compute_LICORS_loglik}
\alias{compute_LICORS_loglik}
\title{Log-likelihood of LICORS model}
\usage{
compute_LICORS_loglik(weight_matrix = NULL, FLC_pdfs = NULL, lambda = 0, penalty = "entropy", 
    q = 2, base = exp(1))
}
\arguments{
  \item{weight_matrix}{\eqn{N \times K} weight matrix}

  \item{FLC_pdfs}{an \eqn{N \times K} matrix containing the
  estimates of all \eqn{K} FLC densities evaluated at all
  \eqn{N} sample FLCs.}

  \item{lambda}{regularization parameter. Default:
  \code{lambda=0} (\code{penalty} and \code{q} will be
  ignored in this case).}

  \item{penalty}{type of penalty: \code{c("entropy",
  "1-Lq", "lognorm")}. Default: \code{"entropy"}}

  \item{base}{logarithm base for the \code{"entropy"}
  penalty. Default: \code{base = 2}.  Any other real number
  is allowed; if \code{base = "nstates"} then it will
  internally assign it \code{base = ncol(weight_matrix)}.}

  \item{q}{exponent for \eqn{L_q} norm.}
}
\description{
  Computes the \emph{average} log-likelihood
  \eqn{\frac{1}{N} \ell(\mathbf{W}; \mathcal{D})} as a
  function of the weight matrix \eqn{\mathbf{W}} and the
  predictive state distributions \eqn{P(X = x \mid S = s_j)
  \approx P(X = x \mid \mathbf{W}_j)} for all \eqn{j = 1,
  \ldots, K}. See References.
}
\keyword{manip}
\keyword{nonparametric}
