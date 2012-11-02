\name{predict_FLC_given_PLC}
\alias{predict_FLC_given_PLC}
\title{Predict FLCs given new PLCs}
\usage{
predict_FLC_given_PLC(state_vector = NULL, FLC_train, FLCs_train_pdf = NULL, 
    type = "mean", weight_matrix_train = NULL, weight_matrix_test = weight_matrix_train)
}
\arguments{
  \item{FLC_train}{the matrix of training FLCs. This matrix
  has \eqn{N} rows, and \eqn{n_f} columns}

  \item{FLCs_train_pdf}{\eqn{N \times K} matrix with
  density estimates of the training FLCs.}

  \item{weight_matrix_train}{\eqn{N \times K} weight matrix
  from the training data}

  \item{weight_matrix_test}{\eqn{\tilde{N} \times K} weight
  matrix from the test data}

  \item{state_vector}{vector of length \eqn{\tilde{N}} with
  entry \eqn{i} being the label \eqn{k = 1, \ldots, K} of
  PLC \eqn{i}}

  \item{type}{type of prediction: \code{'mean'},
  \code{'median'}, \code{'weightedmean'}.}
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
