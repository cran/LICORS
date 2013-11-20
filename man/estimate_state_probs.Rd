\name{estimate_state_probs}
\alias{estimate_state_probs}
\title{Estimate conditional/marginal state probabilities}
\usage{
estimate_state_probs(weight.matrix = NULL, states = NULL, pdfs = list(FLC = NULL, 
    PLC = NULL), num.states = NULL)
}
\arguments{
  \item{weight.matrix}{\eqn{N \times K} weight matrix}

  \item{states}{vector of length \eqn{N} with entry \eqn{i}
  being the label \eqn{k = 1, \ldots, K} of PLC \eqn{i}}

  \item{pdfs}{a list with estimated pdfs for PLC and/or FLC
  evaluated at each PLC, \eqn{i=1, \ldots, N} and/or FLC,
  \eqn{i=1, \ldots, N}}

  \item{num.states}{number of states in total. If
  \code{NULL} (default) then it sets it to
  \code{max(states)} or \code{ncol(weight.matrix)} -
  depending on which one is provided.}
}
\value{
  A vector of length \eqn{K} or a \eqn{N \times K} matrix.
}
\description{
  Estimates \eqn{P(S = s_k; \mathbf{W})}, \eqn{k = 1,
  \ldots, K}, the probability of being in state \eqn{s_k}
  using the weight matrix \eqn{\mathbf{W}}.

  These probabilites can be marginal (\eqn{P(S = s_k;
  \mathbf{W})}) or conditional (\eqn{P(S = s_k \mid
  \ell^{-}, \ell^{+}; \mathbf{W})}), depending on the
  provided information (\code{pdfs$PLC} and
  \code{pdfs$FLC}). \itemize{ \item If both are \code{NULL}
  then \code{estimate_state_probs} returns a vector of
  length \eqn{K} with marginal probabilities.  \item If
  either of them is not \code{NULL} then it returns an
  \eqn{N \times K} matrix, where row \eqn{i} is the
  probability mass function of PLC \eqn{i} being in state
  \eqn{s_k}, \eqn{k = 1, \ldots, K}. }
}
\examples{
WW <- matrix(runif(10000), ncol = 10)
WW <- normalize(WW)
estimate_state_probs(WW)
}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
