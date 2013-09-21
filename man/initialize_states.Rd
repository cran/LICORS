\name{initialize_states}
\alias{initialize_states}
\title{State initialization for iterative algorithms (randomly or variants of kmeans)}
\usage{
initialize_states(nstates = NULL, nsamples = NULL, method = "random", 
    LCs = list(PLC = NULL, FLC = NULL))
}
\arguments{
  \item{nstates}{number of states}

  \item{nsamples}{number of samples.}

  \item{method}{how to choose the labels: either uniformly
  at random from \eqn{\lbrace 1, \ldots, K \rbrace} or
  using K-means on PLCs and FLCs or a combination.
  Default: \code{method = "random"}.  Other options are
  \code{c("KmeansPLC","KmeansFLC","KmeansPLCFLC","KmeansFLCPLC")}}

  \item{LCs}{(optional) a list of \code{PLC} (\eqn{N \times
  n_p} array) and \code{FLC} (\eqn{N \times n_f} array)}
}
\description{
  Initializes the state/cluster assignment either uniformly
  at random from \eqn{K} classes, or using initial
  \emph{kmeans++} (\code{\link{kmeanspp}}) clustering (in
  several variations on PLCs and/or FLCs).
}
\examples{
x1 <- rnorm(1000)
x2 <- rnorm(200, mean = 2)
yy <- c(x1, x2)
ss <- initialize_states(nstates = 2, nsamples = length(yy), method = "KmeansFLC", 
    LCs = list(FLCs = yy))
plot(yy, col = ss, pch = 19)
points(x1, col = "blue")
}
\keyword{datagen}
\keyword{distribution}
\keyword{multivariate}
