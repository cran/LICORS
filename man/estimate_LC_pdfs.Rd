\name{estimate_LC_pdfs}
\alias{estimate_LC_pdfs}
\alias{estimate_LC_pdf_state}
\alias{estimate_LC.pdf.state}
\title{Estimate PLC/FLC distributions for all states}
\usage{
estimate_LC_pdfs(LCs, weight.matrix = NULL, method = c("nonparametric", "normal", 
    "huge"), eval.LCs = NULL)

estimate_LC_pdf_state(state, states = NULL, weights = NULL, LCs = NULL, eval.LCs = NULL, 
    method = c("nonparametric", "normal", "huge"))
}
\arguments{
  \item{LCs}{matrix of PLCs/FLCs. This matrix has \eqn{N}
  rows and \eqn{n_p} or \eqn{n_f} columns (depending on the
  PLC/FLC dimensionality)}

  \item{weight.matrix}{\eqn{N \times K} weight matrix}

  \item{states}{vector of length \eqn{N} with entry \eqn{i}
  being the label \eqn{k = 1, \ldots, K} of PLC \eqn{i}}

  \item{method}{type of estimation: either a (multivariate)
  Normal distribution (\code{"normal"}) or nonparametric
  with a kernel density estimator (\code{method =
  "nonparametric"}). For multivariate distributions (as
  usual for PLCs) only \code{'normal'} should be used due
  to computational efficiency and statistical accuracy.}

  \item{eval.LCs}{on what LCs should the estimate be
  evaluated? If \code{NULL} then densities will be
  evaluated on the training data \code{LCs}}

  \item{state}{integer; which state-conditional density
  should be estimated}

  \item{weights}{weights of the samples. Either a i) length
  \eqn{N} vector with the weights for each observation; ii)
  \eqn{N \times K} matrix, where the \code{state} column of
  that matrix is used as a weight-vector.}
}
\value{
  \code{\link{estimate_LC_pdfs}} returns an \eqn{N \times
  K} matrix.

  \code{\link{estimate_LC.pdf.state}} returns a vector of
  length \eqn{N} with the state-conditional density
  evaluated at \code{eval.LCs}.
}
\description{
  \code{\link{estimate_LC_pdfs}} estimates the PLC and FLC
  distributions for each state \eqn{k = 1, \ldots, K}. It
  iteratively applies \code{\link{estimate_LC.pdf.state}}.

  \code{\link{estimate_LC.pdf.state}} estimates the PLC and
  FLC distributions using weighted maximum likelihood
  (\code{\link[stats]{cov.wt}}) and nonparametric kernel
  density estimation (\code{\link{wKDE}}) for one (!)
  state.
}
\examples{
set.seed(10)
WW <- matrix(runif(10000), ncol = 10)
WW <- normalize(WW)
temp_flcs <- cbind(sort(rnorm(nrow(WW))))
temp_flc_pdfs <- estimate_LC_pdfs(temp_flcs, WW)
matplot(temp_flcs, temp_flc_pdfs, col = 1:ncol(WW), type = "l", xlab = "FLCs", 
    ylab = "pdf", lty = 1)
###################### one state only ###
temp_flcs <- temp_flcs[order(temp_flcs)]
temp_flc_pdf <- estimate_LC_pdf_state(state = 3, LCs = temp_flcs, weights = WW)

plot(temp_flcs, temp_flc_pdf, type = "l", xlab = "FLC", ylab = "pdf")
}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
