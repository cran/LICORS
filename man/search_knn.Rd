\name{search_knn}
\alias{search_knn}
\title{K nearest neighbor (KNN) search}
\usage{
search_knn(data, k = 1, query = NULL, method = "FNN", ...)
}
\arguments{
  \item{data}{an \eqn{N \times d} matrix, where \eqn{N} are
  the samples and \eqn{d} is the dimension of space.  For
  large \eqn{d} knn search can be very slow.}

  \item{k}{number of nearest neighbors (excluding point
  itself). Default: \code{k=1}.}

  \item{query}{(optional) an \eqn{\tilde{N} \times d}
  matrix to find KNN in the training data for. Must have
  the same \eqn{d} as \code{data}; can have lower or larger
  \eqn{\tilde{N}} though. Default: \code{query=NULL}
  meaning that nearest neighbors should be looked for in
  the training data itself.}

  \item{method}{what method should be used: \code{'FNN'},
  \code{'RANN'}, or \code{'yaImpute'}.}

  \item{...}{other parameters passed to the knn functions
  in each package.}
}
\description{
  This is a wrapper for several k nearest neighbors (KNNs)
  algorithms in R. Currently wrapped functions are from the
  \code{FNN}, \code{RANN}, and \code{yaImpute} package.

  It searches for KNN in a \eqn{N \times d} data matrix
  \code{data} where \eqn{N} are the number of samples, and
  \eqn{d} is the dimension of space.

  Either knn search in itself \code{query=NULL} or to query
  new data points wrt to training dataset.
}
\examples{
set.seed(1984)
XX <- matrix(rnorm(40), ncol = 2)
YY <- matrix(runif(length(XX) * 2), ncol = ncol(XX))
knns_of_XX_in_XX <- search_knn(XX, 1)
knns_of_YY_in_XX <- search_knn(XX, 1, query = YY)
plot(rbind(XX, YY), type = "n", xlab = "", ylab = "")
points(XX, pch = 19, cex = 2, xlab = "", ylab = "")
arrows(XX[, 1], XX[, 2], XX[knns_of_XX_in_XX, 1], XX[knns_of_XX_in_XX, 2], lwd = 2)
points(YY, pch = 15, col = 2)
arrows(YY[, 1], YY[, 2], XX[knns_of_YY_in_XX, 1], XX[knns_of_YY_in_XX, 2], col = 2)
legend("left", c("X", "Y"), lty = 1, pch = c(19, 15), cex = c(2, 1), col = c(1, 2))
}
\seealso{
  Packages \pkg{FNN}, \pkg{RANN}, and \pkg{yaImpute} for
  other options (\code{...}).
}
\keyword{classif}
\keyword{cluster}
\keyword{nonparametric}
