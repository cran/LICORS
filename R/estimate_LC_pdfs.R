#' @title Estimate PLC/FLC distributions for all states
#' @aliases estimate_LC_pdf_given_state
#' @description  
#' \code{\link{estimate_LC_pdfs}} estimates the PLC and FLC distributions for 
#' each state \eqn{k = 1, \ldots, K}. It iteratively applies 
#' \code{\link{estimate_LC_pdf_given_state}}.
#' 
#' \code{\link{estimate_LC_pdf_given_state}} estimates the PLC and FLC 
#' distributions using weighted maximum likelihood (\code{\link[stats]{cov.wt}}) 
#' and nonparametric kernel density estimation (\code{\link{wKDE}}) for one
#'  (!) state.
#' 
#' @param LCs matrix of PLCs/FLCs. This matrix has \eqn{N} rows and 
#' \eqn{n_p} or \eqn{n_f} columns (depending on the PLC/FLC dimensionality)
#' @param weight_matrix \eqn{N \times K} weight matrix
#' @param state_vector vector of length \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param method type of estimation: either a (multivariate) Normal distribution 
#' (\code{"normal"}) or nonparametric with a 
#' kernel density estimator (\code{method = "nonparametric"}).  
#' For multivariate distributions (as usual for PLCs) only
#' \code{'normal'} should be used due to computational efficiency and 
#' statistical accuracy.
#' @param eval_LCs on what LCs should the estimate be evaluated? If \code{NULL} then densities
#' will be evaluated on the training data \code{LCs}
#' @return
#' \code{\link{estimate_LC_pdfs}} returns an \eqn{N \times K} matrix.
#' 
#' @keywords nonparametric multivariate distribution
#' @export
#' @examples
#' set.seed(10)
#' WW = matrix(runif(10000), ncol = 10)
#' WW = normalize(WW)
#' temp_flcs = cbind(sort(rnorm(nrow(WW))))
#' temp_flc_pdfs = estimate_LC_pdfs(temp_flcs, WW)
#' matplot(temp_flcs, temp_flc_pdfs, col = 1:ncol(WW), type = "l", 
#'         xlab = "FLCs", ylab = "pdf", lty = 1)
#' 
#' 
estimate_LC_pdfs <- function(LCs, weight_matrix = NULL,
                             state_vector = NULL, method = "normal", 
                             eval_LCs = NULL) {
  
  if (is.null(weight_matrix) & is.null(state_vector)){
    state_vector = rep(1, nrow(LCs))
    weight_matrix = NULL
    nstates <- 1
  } else {  
    if (!is.null(weight_matrix)) {
      nstates <- ncol(weight_matrix)
    } else {
      nstates <- length(unique(state_vector))
    }
  }
  
  if (is.null(eval_LCs)){
    nevals <- nrow(LCs)
  } else {
    nevals <- nrow(eval_LCs)
  }
  pdf_LCs <- matrix(0, ncol = nstates, nrow = nevals)
  for (ii in 1:nstates) {
    pdf_LCs[, ii] <- estimate_LC_pdf_given_state( state_label = ii, 
                                                  state_vector = state_vector,
                                                  weights = weight_matrix[, ii],
                                                  eval_LCs = eval_LCs,
                                                  LCs = LCs, 
                                                  method = method)
  }
  invisible(pdf_LCs)
} 


#' @rdname estimate_LC_pdfs
#' @param state_label integer; which state-conditional density should be 
#' estimated
#' @param weights weights of the samples. Either a i) length \eqn{N} vector with the 
#' weights for each observation; ii) \eqn{N \times K} matrix, where the
#' \code{state_label} column of that matrix is used as a weight-vector.
#' @keywords nonparametric multivariate distribution
#' @return
#' \code{\link{estimate_LC_pdf_given_state}} returns a vector of length \eqn{N} 
#' with the state-conditional density evaluated at \code{eval_LCs}.
#' @export
#' @examples
#' 
#' ######################
#' ### one state only ###
#' ######################
#' temp_flc_pdf = estimate_LC_pdf_given_state(state_label = 3, 
#'                                             LCs = temp_flcs, 
#'                                             weights = WW)
#' plot(temp_flcs, temp_flc_pdf, type = "l", xlab = "FLC", ylab = "pdf")

estimate_LC_pdf_given_state <- function(state_label = NULL, state_vector = NULL, 
                                        weights = NULL, 
                                        LCs = NULL, 
                                        eval_LCs = NULL, 
                                        method = "nonparametric") {
  
  if (is.null(LCs)){
    stop("You have to provide data (LCs) to estimate the densities.")
  }
  
  if (is.null(state_vector) && is.null(weights)) {
    stop("Either provide state assignment or weights.")
  }
  
  if (is.matrix(weights)){
    # necessary below for optimal bandwidth selection
    state_vector <- weight_matrix2states(weights)
    weights <- weights[, state_label]
  }  
  
  # state adaptive bandwidth selection adjust the bw.nrd0() rule to only 
  # observations in the given state
  LCs_in_state <- which(state_vector == state_label)
  if (length(LCs_in_state) < 5){
    LCs_in_state <- which(weights >= quantile(weights, 0.9) )
  }
  optimal_bw <- bw.nrd0(c(LCs)[LCs_in_state])
  
  LCs <- cbind(LCs)
  if (is.null(eval_LCs)) {
    eval_LCs <- LCs
  }
  eval_LCs <- cbind(eval_LCs)
  nsamples <- nrow(eval_LCs)
  kk <- ncol(eval_LCs)
  
  effective_nsamples <- c()
  if (is.null(weights)) {
    effective_nsamples <- sum(state_vector == state_label)
  } else {
    effective_nsamples <- sum(weights)
  }
  
  
  if (kk == 1) {
    # univariate density estimation
    if (method == "nonparametric") {
      if (is.null(weights) || is.na(weights)) {
        LC_pdf_given_state <- wKDE(LCs[LCs_in_state, ], 
                                   eval_points = eval_LCs, 
                                   weights = NULL,
                                   bw = optimal_bw)
      } else {
        LC_pdf_given_state <- wKDE(LCs, eval_points = eval_LCs, 
                                   weights = weights/sum(weights),
                                   bw = optimal_bw)
      }
    } else if (method == "normal") {
      if (is.null(weights) || is.na(weights)) {
        mu <- mean(LCs[state_vector == state_label, ])
        sigma <- sd(LCs[state_vector == state_label, ])
      } else {
        temp <- cov.wt(LCs, wt = weights, method = "ML")
        mu <- temp$center
        sigma <- sqrt(temp$cov)
        rm(temp)
      }
      
      if (floor(effective_nsamples) < kk) {
        sigma <- sd(LCs)
      }
      
      LC_pdf_given_state <- dnorm(c(eval_LCs), mean = mu, sd = sigma)
    }
  } else {
    # Multivariate density estimation
    if (method == "nonparametric") {
      stop('Multivariate nonparametrc not implemented yet (too CPU intensive 
           for high dimensions).')
      
      if (is.null(weights) || is.na(weights)) {
        LC_pdf_given_state <- mv_wKDE(LCs[state_vector == state_label, ], 
                                      eval_points = eval_LCs)
      } else {
        LC_pdf_given_state <- mv_wKDE(LCs, eval_points = eval_LCs, 
                                      weights = weights/sum(weights))
      }
    } else if (method == "normal") {
      
      if (floor(effective_nsamples) < kk) {
        Sigma_cov <- cov(LCs)
        mu_vec <- apply(LCs, 2, median)
      } else {
        if (is.null(weights) || is.na(weights)) {
          mu_vec <- colMeans(LCs[state_vector == state_label, ])
          Sigma_cov <- cov(LCs[state_vector == state_label, ])
        } else {
          temp <- cov.wt(LCs, wt = weights, method = "ML")
          mu_vec <- temp$center
          Sigma_cov <- temp$cov
          rm(temp)
        }
      }
      Sigma_cov <- Sigma_cov + 1/100 * mean(diag(Sigma_cov)) * diag(1, ncol(Sigma_cov))
      LC_pdf_given_state <- dmvnorm(eval_LCs, mean = mu_vec, sigma = Sigma_cov)
    } else if (method == "huge") {
      glasso.est <- huge.select(huge(LCs[state_vector == state_label, ], 
                                     cov.output = TRUE, method = "glasso", 
                                     verbose = FALSE), 
                                verbose = FALSE)
      LC_pdf_given_state <- dmvnorm(eval_LCs, 
                                    mean = apply(LCs[state_vector == state_label, ], 2, median), 
                                    sigma = as.matrix(glasso.est$opt.cov))
    }
    
  }
  LC_pdf_given_state[is.na(LC_pdf_given_state)] <- .Machine$double.eps^0.25
  invisible(LC_pdf_given_state)
}


