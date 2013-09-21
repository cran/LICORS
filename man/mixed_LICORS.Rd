\name{mixed_LICORS}
\alias{mixed_LICORS}
\title{Mixed LICORS: An EM-like Algorithm for Predictive State Space Estimation}
\usage{
mixed_LICORS(LCs = list(PLC = PLCs, FLC = FLCs, dim = list(original = NULL, 
    truncated = NULL)), nstates_start = NULL, initialization = NULL, max_iter = 500, 
    method = list(PLC = "normal", FLC = "nonparametric"), alpha = 0.01, 
    seed = NULL, verbose = TRUE, CV_train_ratio = 0.75, CV_split_random = FALSE, 
    sparsity = "stochastic", lambda = 0, loss = function(x, xhat) mean((x - 
        xhat)^2))
}
\arguments{
  \item{LCs}{list of PLCs and FLCs matrices (see output of
  \code{\link{data2LCs}} for details and formatting).}

  \item{nstates_start}{number of states to start the EM
  algorithm}

  \item{initialization}{a a) character string, b) vector,
  or c) matrix. a) results \code{nstates_start} many states
  initialized by passing the character string as
  \code{method} argument of
  \code{\link{initialize_states}}; if b) the vector will be
  taken as initial state labels; if c) the matrix will be
  taken as initial weights. Note that for both b) and c)
  \code{nstates_start} will be ignored. \eqn{k = 1, \ldots,
  K} of PLC \eqn{i}}

  \item{max_iter}{maximum number of iterations in the EM}

  \item{verbose}{logical; if \code{TRUE} it prints output
  in the console as the EM is running}

  \item{sparsity}{what type of sparsity (currently not
  implemented)}

  \item{lambda}{penalization parameter; larger lambda gives
  sparser weights}

  \item{loss}{an R function specifying the loss for
  cross-validation (CV). Default: mean squared error (MSE),
  i.e. \code{loss = function(x, xhat) mean((x-xhat)^2)}}

  \item{alpha}{significance level to stop testing. Default:
  \code{alpha=0.01}}

  \item{seed}{set seed for reproducibility. Default:
  \code{NULL}. If \code{NULL} it sets a random seed and
  then returns this seed in the output.}

  \item{CV_train_ratio}{how much of the data should be
  training data. Default: \code{0.75}, i.e. \eqn{75\%} of
  data is for training}

  \item{CV_split_random}{indicator for splitting data
  randomly in training and test data (\code{TRUE}) or to
  use first part (in time) as training, rest as test
  (\code{FALSE}; default).}

  \item{method}{a list of length \eqn{2} with arguments
  \code{PLC} and \code{FLC} for the method of density
  estimation in each (either \code{"normal"} or
  \code{"nonparametric"}).}
}
\value{
  An object of class \code{"LICORS"}.
}
\description{
  \code{mixed_LICORS} is the core function of this package
  as it estimates the ``parameters'' in the model for the
  spatio-temporal process. \deqn{ P(X_1, \ldots,
  X_{\tilde{N}}) \propto \prod_{i=1}^{N} P(X_i \mid
  \ell^{-}_i) = \prod_{i=1}^{N} P(X_i \mid
  \epsilon(\ell^{-}_i)) . }
}
\examples{
\dontrun{
data(contCA00)

LC_geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 2, FLC = 0), 
    shape = "cone")
bb <- data2LCs(t(contCA00$observed), LC_coords = LC_geom$coords)

mm <- mixed_LICORS(bb, nstates_start = 10, init = "KmeansPLC", max_iter = 20)
plot(mm)
ff_new <- estimate_LC_pdfs(bb$FLC, weight_matrix = mm$conditional_state_probs, 
    method = "nonparametric")
matplot(bb$FLC, ff_new, pch = ".", cex = 2)
}
}
\seealso{
  \code{\link{plot.mixed_LICORS}},
  \code{\link{summary.mixed_LICORS}}
}
\keyword{cluster}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
