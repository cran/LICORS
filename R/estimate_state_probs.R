#' @title Estimate conditional/marginal state probabilities
#'
#' @description 
#' Estimates \eqn{P(S = s_k; \mathbf{W})}, \eqn{k = 1, \ldots, K}, the probability of being 
#' in state \eqn{s_k} using the weight matrix \eqn{\mathbf{W}}.
#' 
#' These probabilites can be marginal (\eqn{P(S = s_k; \mathbf{W})}) or 
#' conditional (\eqn{P(S = s_k \mid \ell^{-}, \ell^{+}; \mathbf{W})}), depending on the 
#' provided information (\code{PLC_pdfs} and \code{FLC_pdfs}).  
#' \itemize{
#'   \item If both are \code{NULL} then \code{estimate_state_probs} returns a 
#'         vector of length \eqn{K} with marginal probabilities.
#'   \item If either of them is not \code{NULL} then it returns an 
#'         \eqn{N \times K} matrix, where row \eqn{i} is the 
#'         probability mass function of PLC \eqn{i} being in state 
#'         \eqn{s_k}, \eqn{k = 1, \ldots, K}.
#' }
#'  
#' @param weight_matrix \eqn{N \times K} weight matrix
#' @param state_vector vector of length \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param PLC_pdfs estimated PLC pdf evaluated at each PLC, \eqn{i=1, \ldots, N}
#' @param FLC_pdfs estimated FLC pdf evaluated at each FLC, \eqn{i=1, \ldots, N}
#' @param nstates_total number of states in total. If \code{NULL} (default) then
#' it sets it to \code{max(state_vector)} or \code{ncol(weight_matrix)} - 
#' depending on which one is provided.
#' @return
#' A vector of length \eqn{K} or a \eqn{N \times K} matrix.
#' @keywords nonparametric multivariate distribution
#' @export
#' @examples
#' WW = matrix(runif(10000), ncol = 20)
#' WW = normalize(WW)
#' estimate_state_probs(WW)

estimate_state_probs <- function(weight_matrix = NULL, 
                                 state_vector= NULL, 
                                 PLC_pdfs = NULL, 
                                 FLC_pdfs = NULL, 
                                 nstates_total = NULL) {
  
  if (is.null(state_vector) & is.null(weight_matrix)){
    stop("You must provide either a vector of state 
          labels or the weight matrix.")
  }
  

  if (is.null(PLC_pdfs) & is.null(FLC_pdfs)){
    ##############################
    ### marginal probabilities ###
    ##############################
    
    if (is.null(nstates_total)){
      if (is.null(state_vector)){
        nstates_total <- ncol(weight_matrix)
      } else {
        nstates_total <- max(state_vector)
      }
      
    }
    
    if (is.null(weight_matrix)) {
      state_probs = c((table(c(state_vector, 1:nstates_total)) - 1)) / length(state_vector)
    } else {
      state_probs = colMeans(weight_matrix)
    }
    
  } else { 
    #################################
    ### conditional probabilities ###
    #################################
    marginal_state_probs <- estimate_state_probs(weight_matrix = weight_matrix, 
                                                 state_vector = state_vector,
                                                 nstates_total = nstates_total)
    
    if (is.null(nstates_total) && is.null(PLC_pdfs) && is.null(FLC_pdfs)) {
      stop("you must provide either one of those")
    }
    
    if (is.null(nstates_total)) {
      nstates_total <- ncol(PLC_pdfs)
    }
    # initialize pi_theta
    if (is.null(PLC_pdfs)) {
      state_probs_given_LCs <- matrix(0, ncol = nstates_total, nrow = length(state_vector))
    } else {
      state_probs_given_LCs <- matrix(0, ncol = nstates_total, nrow = nrow(PLC_pdfs))
    }
    
    # if PLCs are not provided, initialize weights with 0/1 assignment based on
    # state
    if (is.null(PLC_pdfs)) {
      state_probs_given_LCs[cbind(1:nrow(state_probs_given_LCs), state_vector)] = 1
    } else {
      state_probs_given_LCs <- sweep(PLC_pdfs, 2, marginal_state_probs, "*")
    } 
    # condition on FLCs if provided
    if (!is.null(FLC_pdfs)) {
      state_probs_given_LCs <- state_probs_given_LCs * FLC_pdfs
    }
  }
  
  if (is.null(PLC_pdfs) & is.null(FLC_pdfs)){
    return(state_probs)
  } else {
    invisible(normalize(state_probs_given_LCs))
  }
} 