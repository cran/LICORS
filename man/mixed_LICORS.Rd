\name{mixed_LICORS}
\alias{mixed_LICORS}
\title{Mixed LICORS: An EM-like Algorithm for Predictive State Space Estimation}
\usage{
mixed_LICORS(LCs = list(PLC = NULL, FLC = NULL, dim = list(original = NULL, 
    truncated = NULL)), num.states.init = NULL, initialization = NULL, 
    control = list(max.iter = 500, alpha = 0.01, trace = 0, lambda = 0, 
        sparsity = "stochastic", CV.split.random = FALSE, CV.train.ratio = 0.75, 
        seed = NULL, loss = function(x, xhat) mean((x - xhat)^2), 
        estimation.method = list(PLC = "normal", FLC = "nonparametric")))
}
\arguments{
  \item{LCs}{list of PLCs and FLCs matrices (see output of
  \code{\link{data2LCs}} for details and formatting).}

  \item{num.states.init}{number of states to start the EM
  algorithm}

  \item{initialization}{a a) character string, b) vector,
  or c) matrix. a) results \code{num.states.init} many
  states initialized by passing the character string as
  \code{method} argument of
  \code{\link{initialize_states}}; if b) the vector will be
  taken as initial state labels; if c) the matrix will be
  taken as initial weights. Note that for both b) and c)
  \code{num.states.init} will be ignored.  \eqn{k = 1,
  \ldots, K} of PLC \eqn{i}}

  \item{control}{a list of control settings for the EM
  algorithm. See \code{\link{complete_LICORS_control}} for
  details.}
}
\value{
  An object of class \code{"LICORS"}.
}
\description{
  \code{mixed_LICORS} is the core function of this package
  as it estimates the ``parameters'' in the model for the
  spatio-temporal process.  \deqn{ P(X_1, \ldots,
  X_{\tilde{N}}) \propto \prod_{i=1}^{N} P(X_i \mid
  \ell^{-}_i) = \prod_{i=1}^{N} P(X_i \mid
  \epsilon(\ell^{-}_i)) .  }
}
\examples{
\dontrun{
data(contCA00)

LC_geom <- setup_LC_geometry(speed = 1, horizon = list(PLC = 2, FLC = 0), 
    shape = "cone")
bb <- data2LCs(t(contCA00$observed), LC.coordinates = LC_geom$coordinates)

mm <- mixed_LICORS(bb, num.states.init = 15, init = "KmeansPLC", 
    control = list(max.iter = 50, lambda = 0.001))
plot(mm)
ff_new <- estimate_LC_pdfs(bb$FLC, weight.matrix = mm$conditional_state_probs, 
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
