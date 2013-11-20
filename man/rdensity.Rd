\name{rdensity}
\alias{rdensity}
\title{Generate random sample from density() or wKDE}
\usage{
rdensity(n = 100, data = NULL, fhat = NULL, bw = fhat$bw, weights = NULL, 
    kernel = "Gaussian")
}
\arguments{
  \item{n}{number of samples}

  \item{fhat}{an object of class
  \code{'\link[stats]{density}'} for bandwidth selection
  (if \code{bw} is not explicitly provided as argument)}

  \item{weights}{vector of weights. Same length as
  \code{data}. Default \code{weights=NULL} - in this case
  equal weights for each point}

  \item{data}{underlying sample of \code{fhat}}

  \item{kernel}{kernel choice for \code{fhat}. Default:
  \code{kernel='Gaussian'}. See
  \code{\link[stats]{density}} for other options.}

  \item{bw}{choice of bandwidth. Default:
  \code{bw=fhat$bw}. Again see \code{\link[stats]{density}}
  for other options.}
}
\description{
  This function draws random samples given \code{data} and
  a \code{\link[stats]{density}} estimate (or just
  providing the correct bandwidth \code{bw}).
}
\examples{
set.seed(1923)
xx <- c(rnorm(100, mean = 2), runif(100))
aa <- density(xx)
plot(aa)
xx_sample <- rdensity(n = 1000, fhat = aa, data = xx)
lines(density(xx_sample), col = 2)
}
\keyword{#}
\keyword{distribution}
\keyword{nonparametric}
