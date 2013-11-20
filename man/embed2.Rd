\name{embed2}
\alias{embed2}
\title{Improved embed() function}
\usage{
embed2(x, max.lag = 1, na.omit = FALSE)
}
\arguments{
  \item{x}{a numeric vector, matrix, or time series.}

  \item{max.lag}{a scalar representing the embedding
  dimension in past or future. Note that contrary to
  'dimension = 1' in \code{\link[stats]{embed}}, here
  'max.lag = 1' will return a 2 column matrix (0 lag, 1
  lag), and not just a 1 column matrix.  Similarly, for
  negative shift; e.g., 'max.lag = -2' returns 3 column
  matrix with (0 lag, -1 lag, -2 lag).}

  \item{na.omit}{logical; if TRUE, it removes NA values
  automatically from embedded matrix.}
}
\description{
  Improved version of the \code{\link[stats]{embed}}
  function in the \code{stats} package. First it allows
  embeddings in past and future observation space (backward
  and forward shift). Secondly, it adds 'NA' to the
  beginning (or end) of the embedding matrix, depending on
  the dimension of the embedding. Optionally, they can be
  removed.
}
\examples{
data(nottem)
aa <- embed2(nottem, 12)
}
\seealso{
  \code{\link[stats]{embed}}
}
