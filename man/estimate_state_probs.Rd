\name{estimate_state_probs}
\alias{estimate_state_probs}
\title{Estimate conditional/marginal state probabilities}
\usage{
estimate_state_probs(weight_matrix = NULL, state_vector = NULL, PLC_pdfs = NULL, 
    FLC_pdfs = NULL, nstates_total = NULL)
}
\arguments{
  \item{weight_matrix}{\eqn{N \times K} weight matrix}

  \item{state_vector}{vector of length \eqn{N} with entry
  \eqn{i} being the label \eqn{k = 1, \ldots, K} of PLC
  \eqn{i}}

  \item{PLC_pdfs}{estimated PLC pdf evaluated at each PLC,
  \eqn{i=1, \ldots, N}}

  \item{FLC_pdfs}{estimated FLC pdf evaluated at each FLC,
  \eqn{i=1, \ldots, N}}

  \item{nstates_total}{number of states in total. If
  \code{NULL} (default) then it sets it to
  \code{max(state_vector)} or \code{ncol(weight_matrix)} -
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
  provided information (\code{PLC_pdfs} and
  \code{FLC_pdfs}). \itemize{ \item If both are \code{NULL}
  then \code{estimate_state_probs} returns a vector of
  length \eqn{K} with marginal probabilities.  \item If
  either of them is not \code{NULL} then it returns an
  \eqn{N \times K} matrix, where row \eqn{i} is the
  probability mass function of PLC \eqn{i} being in state
  \eqn{s_k}, \eqn{k = 1, \ldots, K}. }
}
\examples{
WW <- matrix(runif(10000), ncol = 20)
WW <- normalize(WW)
estimate_state_probs(WW)
}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
