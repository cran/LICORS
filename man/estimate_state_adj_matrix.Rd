\name{estimate_state_adj_matrix}
\alias{estimate_state_adj_matrix}
\title{Estimate adjacency matrix for equivalent FLC distributions based on states}
\usage{
estimate_state_adj_matrix(state_vector = NULL, FLCs = NULL, FLC_pdfs = NULL, 
    alpha = NULL, distance = function(f, g) return(mean(abs(f - g))))
}
\arguments{
  \item{state_vector}{vector of length \eqn{N} with entry
  \eqn{i} being the label \eqn{k = 1, \ldots, K} of PLC
  \eqn{i}}

  \item{FLCs}{\eqn{N \times n_f} matrix of FLCs (only
  necessary if \code{distance= "KS"})}

  \item{FLC_pdfs}{\eqn{N \times K} matrix of all \eqn{K}
  state-conditional FLC densities evaluated at each FLC
  \eqn{\ell^{+}_i}, \eqn{i=1, \ldots, N} (only necessary if
  \code{distance = function(f, g) return(...)}).}

  \item{alpha}{significance level for testing. Default:
  \code{alpha=NULL} (this will return a p-value matrix if
  \code{method == "KS"})}

  \item{distance}{either a Kolmogorov-Smirnov test
  (\code{distance = "KS"}) or a function metric (e.g.
  \eqn{L_q} distance). For a distance function,
  \code{distance} requires as input a function of \eqn{f}
  and \eqn{g} that returns one value.

  Default: \code{distance = function(f, g)
  return(mean(abs(f-g)))} \eqn{\rightarrow} \eqn{L_1}
  distance.}
}
\value{
  A \eqn{K \times K} adjacency matrix with a trimmed
  version of exp(-distance) or p-values.  If
  \code{alpha!=NULL} then it returns the thresholded
  \eqn{0/1} matrix.  However, here \eqn{1} stands for
  equivalent, i.e. not rejecting. The matrix is obtained by
  checking for \code{pval>alpha} (rather than the usual
  \code{pval<alpha}).
}
\description{
  This function estimates the adjacency matrix
  \eqn{\mathbf{A}} of all pairwise equivalent FLC
  distributions given the states \eqn{s_1, \ldots, s_K}.
  See Details below.
}
\section{Details and user-defined distance function}{
  The \eqn{(i,j)}th element of the adjacency matrix is
  defined as \deqn{ \mathbf{A}_{ij} = distance(P(X \mid
  s_i), P(X \mid s_j)) = distance(f, g), } where
  \code{distance} is either \describe{ \item{a metric}{in
  the function space of pdfs \eqn{f} and \eqn{g}, or}
  \item{a two sample test}{for \eqn{H_0: f=g}, e.g. a
  Kolmogorov-Smirnov test (\code{distance="KS"}).} }

  Again we use a functional programming approach and allow
  the user to specify any valid distance/similarity
  function \code{distance = function(f, g) return(...)}.

  If \code{distance="KS"} the adjacency matrix contains
  p-values of a Kolmogorov-Smirnov test or the thresholded
  versions (if \code{alpha!=NULL}) - see Return for
  details.

  Otherwise \code{distance} is an R function that takes as
  an input two vectors \code{f} and \code{g} (e.g. the
  \code{\link{wKDE}} estimates for two states), and returns
  a non-negative, real number to estimate their distance.
  Default is the \eqn{L_1} distance \code{distance =
  function(f, g) return(mean(abs(f-g)))}.
}
\examples{
WW <- matrix(runif(10000), ncol = 10)
WW <- normalize(WW)
temp_flcs <- cbind(rnorm(nrow(WW)))
temp_flc_pdfs <- estimate_LC_pdfs(temp_flcs, WW)
AA_ks <- estimate_state_adj_matrix(state_vector = weight_matrix2states(WW), 
    FLCs = temp_flcs, distance = "KS")
AA_L1 <- estimate_state_adj_matrix(FLC_pdfs = temp_flc_pdfs)

par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
image2(AA_ks, zlim = c(0, 1), legend = FALSE, main = "Kolmogorov-Smirnov")
image2(AA_L1, legend = FALSE, main = "L1 distance")
}
\keyword{distribution}
\keyword{manip}
\keyword{multivariate}
\keyword{nonparametric}
