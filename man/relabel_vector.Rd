\name{relabel_vector}
\alias{relabel_vector}
\title{Relabels a vector to consecutive labels}
\usage{
relabel_vector(vec, order = FALSE)
}
\arguments{
  \item{vec}{vector with labels}

  \item{order}{logical; if \code{TRUE} then new state
  labels are assigned by decreasing number of points in
  that state. That is, state ``1'' has the most points in
  the state, followed by state ``2'' etc.}
}
\description{
  This function relabels a vector to have consecutive - no
  missing in between - labels.  Labels always start at
  \eqn{1} and increase by one.

  For example, \code{c(2, 2, 5)} gets relabeled to
  \code{c(1, 1, 2)}.
}
\examples{
TempVec <- c(10, 2, 1, 2, 2, 2, 10)
print(relabel_vector(TempVec))

print(relabel_vector(c(2, 2, 5)))
}
\keyword{arith}
\keyword{array}
\keyword{list}
\keyword{manip}
