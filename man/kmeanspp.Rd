\name{kmeanspp}
\alias{kmeanspp}
\title{Kmeans++}
\usage{
kmeanspp(data, k = 2, start = "random", iter.max = 100, nstart = 10, ...)
}
\arguments{
  \item{data}{an \eqn{N \times d} matrix, where \eqn{N} are
  the samples and \eqn{d} is the dimension of space.}

  \item{k}{number of clusters.}

  \item{start}{first cluster center to start with}

  \item{iter.max}{the maximum number of iterations allowed}

  \item{nstart}{how many random sets should be chosen?}

  \item{...}{additional arguments passed to
  \code{\link[stats]{kmeans}}}
}
\description{
  \emph{kmeans++} clustering (see References) using R's
  built-in function \code{\link[stats]{kmeans}}.
}
\examples{
set.seed(1984)
nn <- 100
XX <- matrix(rnorm(nn), ncol = 2)
YY <- matrix(runif(length(XX) * 2, -1, 1), ncol = ncol(XX))
ZZ <- rbind(XX, YY)

cluster_ZZ <- kmeanspp(ZZ, k = 5, start = "random")

plot(ZZ, col = cluster_ZZ$cluster + 1, pch = 19)
}
\references{
  Arthur, D. and S. Vassilvitskii (2007). ``k-means++: The
  advantages of careful seeding.'' In H. Gabow (Ed.),
  Proceedings of the 18th Annual ACM-SIAM Symposium on
  Discrete Algorithms [SODA07], Philadelphia, pp.
  1027-1035. Society for Industrial and Applied
  Mathematics.
}
\seealso{
  \code{\link[stats]{kmeans}}
}
\keyword{cluster}
\keyword{multivariate}
