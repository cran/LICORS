\name{mixed_LICORS-utils}
\alias{complete_LICORS_control}
\alias{mixed_LICORS-utils}
\alias{plot.mixed_LICORS}
\alias{predict.mixed_LICORS}
\alias{summary.mixed_LICORS}
\title{Utilities for ``LICORS'' class}
\usage{
\method{plot}{mixed_LICORS}(x, type = "both", cex.axis = 1.5, cex.lab = 1.5, 
    cex.main = 2, line = 1.5, ...)

\method{summary}{mixed_LICORS}(object, ...)

\method{predict}{mixed_LICORS}(object, new.LCs = list(PLC = NULL), ...)

complete_LICORS_control(control = list(alpha = 0.01, CV.split.random = FALSE, 
    CV.train.ratio = 0.75, lambda = 0, max.iter = 500, seed = NULL, 
    sparsity = "stochastic", trace = 0, loss = function(x, xhat) mean((x - 
        xhat)^2), estimation.method = list(PLC = "normal", FLC = "nonparametric")))
}
\arguments{
  \item{x}{object of class \code{"mixed_LICORS"}}

  \item{type}{should only \code{"training"}, \code{"test"},
  or \code{"both"} be plotted. Default: \code{"both"}.}

  \item{cex.axis}{The magnification to be used for axis
  annotation relative to the current setting of
  \code{cex}.}

  \item{cex.lab}{The magnification to be used for x and y
  labels relative to the current setting of \code{cex}.}

  \item{cex.main}{The magnification to be used for main
  titles relative to the current setting of \code{cex}.}

  \item{line}{on which margin line should the labels be
  ploted, starting at 0 counting objectwards (see also
  \code{\link[graphics]{mtext}}).}

  \item{...}{optional arguments passed to \code{plot},
  \code{summary}, or \code{predict}}

  \item{object}{object of class \code{"mixed_LICORS"}}

  \item{new.LCs}{a list with PLC configurations to predict
  FLCs given these PLCs}

  \item{control}{a list of controls for
  \code{"mixed_LICORS"}.}
}
\description{
  The \code{"mixed_LICORS"} class is the objectput from the
  \code{\link{mixed_LICORS}} estimator.

  \code{plot.mixed_LICORS} gives a visual summary of the
  estimates such as marginal state probabilities,
  conditional state probabilities (= weight matrix),
  predictive state densities, trace plots for
  log-likelihood/loss/penalty.

  \code{summary.mixed_LICORS} prints object a summary of
  the estimated LICORS model.

  \code{predict.mixed_LICORS} predicts FLCs based on PLCs
  given a fitted mixed LICORS model. This can be done on an
  iterative basis, or for a selection of future PLCs.

  \code{complete_LICORS_control} completes the controls for
  the mixed LICORS estimator. Entries of the list are:

  'loss' an R function specifying the loss for
  cross-validation (CV). Default: mean squared error (MSE),
  i.e. \code{loss = function(x, xhat) mean((x-xhat)^2)}

  'method' a list of length \eqn{2} with arguments
  \code{PLC} and \code{FLC} for the method of density
  estimation in each (either \code{"normal"} or
  \code{"nonparametric"}).

  'max.iter' maximum number of iterations in the EM

  'trace' if > 0 it prints output in the console as the EM
  is running

  'sparsity' what type of sparsity (currently not
  implemented)

  'lambda' penalization parameter; larger lambda gives
  sparser weights

  'alpha' significance level to stop testing. Default:
  \code{alpha = 0.01}

  'seed' set seed for reproducibility. Default:
  \code{NULL}. If \code{NULL} it sets a random seed and
  then returns this seed in the output.

  'CV.train.ratio' how much of the data should be training
  data. Default: \code{0.75}, i.e., \eqn{75\%} of data is
  for training

  'CV.split.random' logical; if \code{TRUE} training and
  test data are split randomly; if \code{FALSE} (default)
  it uses the first part (in time) as training, rest as
  test.

  'estimation' a list of length \eqn{2} with arguments
  \code{PLC} and \code{FLC} for the method of density
  estimation in each (either \code{"normal"} or
  \code{"nonparametric"}).
}
\examples{
# see examples of LICORS-package see examples in LICORS-package see examples in
# LICORS-package see examples in LICORS-package
}
\keyword{hplot}
\keyword{model}
\keyword{nonparametric}
