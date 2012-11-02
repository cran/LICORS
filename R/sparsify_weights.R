#' @title Sparsify weights
#'
#' @description 
#' This function makes weights of a mixture model more sparse using gradient 
#' based penalty methods.
#' 
#' @param weight_matrix_proposed \eqn{N \times K} weight matrix
#' @param weight_matrix_current \eqn{N \times K} weight matrix
#' @param penalty type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param lambda penalization parameter: larger \code{lambda} gives sparser 
#' mixture weights
#' @keywords manip array
#' @export
#' @seealso \code{\link{compute_mixture_penalty}}, \code{\link{mixed_LICORS}}
#' @examples
#' WW = matrix(c(rexp(10, 1/10), runif(10)), ncol=5, byrow = FALSE)
#' WW = normalize(WW)
#' WW_sparse = sparsify_weights(WW, lambda = 0.1)
#' WW_more_sparse = sparsify_weights(WW, lambda = 0.5)
#' compute_mixture_penalty(WW)
#' compute_mixture_penalty(WW_sparse)
#' compute_mixture_penalty(WW_more_sparse)

sparsify_weights <- function(weight_matrix_proposed, 
                             weight_matrix_current = NULL, 
                             penalty = "entropy", 
                             lambda = 0) {
  
  if (lambda == 0) {
    return(weight_matrix_proposed)
  } else {
    kk = ncol(weight_matrix_proposed)
    if (is.null(weight_matrix_current)) {
      penalty_gradient = -(1 + log2(weight_matrix_proposed + 10^(-3)))/log2(kk)
    } else {
      penalty_gradient = -(1 + log2(weight_matrix_current + 10^(-3)))/log2(kk)
    }
    # lambda = lambda + runif(1, -lambda/10, lambda/10)
    
    # without hessian weight_matrix_new = weight_matrix_proposed + lambda *
    # penalty_gradient multiply by Hessian
    weight_matrix_new = weight_matrix_proposed * (1 - log2(kk) * lambda * penalty_gradient)
    weight_matrix_new[weight_matrix_new < 0] = 0
    all_zeros = apply(weight_matrix_new, 1, function(u) all(u == 0))
    
    if (sum(all_zeros) > 1) {
      weight_matrix_new[all_zeros, apply(weight_matrix_proposed[all_zeros, ], 1, which.max)] = 1
    } else {
      weight_matrix_new[all_zeros, which.max(weight_matrix_proposed[all_zeros, ])] = 1
    }
    return(normalize(weight_matrix_new))
  }
} 