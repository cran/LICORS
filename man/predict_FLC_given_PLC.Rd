\name{predict_FLC_given_PLC}
\alias{predict_FLC_given_PLC}
\title{Predict FLCs given new PLCs}
\usage{
predict_FLC_given_PLC(train = list(data = list(FLC = NULL, PLC = NULL), 
    weight.matrix = NULL, pdfs = list(FLC = NULL, PLC = NULL)), test = list(PLC = NULL, 
    weight.matrix = NULL), type = c("weighted.mean", "mean", "median", "mode"), 
    method = list(FLC = "nonparametric", PLC = "normal"))
}
\arguments{
  \item{train}{a list of training examples with LC
  observations (a list of \code{PLC} and \code{FLC}),
  \code{weight.matrix}, and \code{pdfs}}

  \item{test}{a list of test examples with PLC observations
  and/or the \code{weight.matrix} associated with the PLC
  observations.}

  \item{method}{estimation method for estimating PLC and
  FLC distributions}

  \item{type}{prediction: \code{'mean'}, \code{'median'},
  \code{'weightedmean'}, or \code{'mode'}.}
}
\value{
  \eqn{N \times K} matrix
}
\description{
  This function predicts FLCs given new PLCs based on the
  estimated \eqn{\epsilon} mappings and estimated
  conditional distributions.
}
\keyword{methods}
