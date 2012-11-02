\name{wKDE}
\alias{mv_wKDE}
\alias{wKDE}
\title{Weighted kernel density estimator (wKDE)}
\usage{
wKDE(x, eval_points = x, weights = NULL, kernel = "gaussian", bw = "nrd0")

mv_wKDE(x, eval_points = x, weights = NULL, kernel = "gaussian")
}
\arguments{
  \item{x}{data vector}

  \item{eval_points}{points where the density should be
  evaluated. Default: \code{eval_points = x}.}

  \item{weights}{vector of weights. Same length as
  \code{x}. Default: \code{weights=NULL} - equal weight for
  each sample.}

  \item{kernel}{type of kernel. Default:
  \code{kernel='Gaussian'}. See
  \code{\link[stats]{density}} and
  \code{\link[locfit]{locfit.raw}} for additional options.}

  \item{bw}{bandwidth. Either a character string indicating
  the method to use or a real number. Default:
  \code{bw="nrd0"}. Again see \code{\link[stats]{density}}
  for other options.}
}
\value{
  A vector of length \code{length(eval_points)} (or
  \code{nrow(eval_points)}) with the probabilities of each
  point given the nonparametric fit on \code{x}.
}
\description{
  \code{wKDE} gives a (weighted) kernel density estimate
  (KDE) for univariate data.

  If weights are not provided, all samples count equally.
  It evaluates on new data point by interpolation (using
  \code{\link[stats]{approx}}).

  \code{mv_KDE} uses the \code{\link[locfit]{locfit.raw}}
  function in the \pkg{locfit} package to estimate KDEs for
  multivariate data. Note: Use this only for small
  dimensions, very slow otherwise.
}
\examples{
### Univariate example ###
xx <- sort(c(rnorm(100, mean = 1), runif(100)))
plot(xx, wKDE(xx), type = "l")
yy <- sort(runif(50, -1, 4) - 1)
lines(yy, wKDE(xx, yy), col = 2)
### Multivariate example ###
XX <- matrix(rnorm(100), ncol = 2)
YY <- matrix(runif(40), ncol = 2)
mm <- mv_wKDE(XX)

plot(mm)
points(mv_wKDE(XX, YY), col = 2, ylab = "")
}
\keyword{distribution}
\keyword{smooth}
