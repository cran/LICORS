#' @title Penalty of mixture weights
#'
#' @description 
#' Computes the penalty \eqn{\Omega(\mathbf{W})} of the weight matrix (or vector) for a mixture model.
#'
#' @param weight_matrix \eqn{N \times K} weight matrix
#' @param type type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param q exponent for \eqn{L_q} norm.
#' @param row_average logical; if \code{TRUE} (default) then an average penalty
#' over all rows will be returned (one single number); if \code{FALSE} a vector 
#' of length \eqn{N} will be returned.
#' @param base logarithm base for the \code{"entropy"} penalty. 
#' Default: \code{base = 2}.  Any other real number is allowed; 
#' if \code{base = "nstates"} then it will internally assign it 
#' \code{base = ncol(weight_matrix)}.
#' @keywords manip array
#' @export
#' @seealso \code{\link{compute_LICORS_loglik}} \code{\link{compute_NEC}}
#' @examples
#' WW = matrix(c(rexp(10, 1/10), runif(10)), ncol=2, byrow = FALSE)
#' WW = normalize(WW)
#' compute_mixture_penalty(WW, row_average = FALSE)
#' compute_mixture_penalty(WW, row_average = TRUE) # default: average penalty

compute_mixture_penalty <- function(weight_matrix, type = "entropy", q = 2, 
                                    row_average = TRUE, base = 2) {
  qq <- q
  kk <- ncol(weight_matrix)
  
  if (base == "nstates"){
    base <- kk
  }
  
  nsamples <- nrow(weight_matrix)
  if (type == "1-Lq") {
    min_norm <- (1/kk^(qq - 1))^(1/qq)
    # f(x) = a x + b; b = -a, and a = 1 / (min_norm(w) - 1))
    aa <- 1/(min_norm - 1)
    penalty <- aa * ((rowSums(weight_matrix^qq))^(1/qq) - 1)
  } else if (type == "entropy") {
    temp <- - weight_matrix * log(weight_matrix, base)
    temp[is.na(temp)] <- 0
    penalty <- rowSums(temp)
  } else if (type == "lognorm") {
    penalty <- -log(rowSums(weight_matrix^qq)^(1/qq))
  } else if (type == "MDL") {
    penalty <- rep(log(kk, base), nsamples)
  } else {
    stop(paste("Penalty of type '", type, "' is not implemented.", sep = ""))
  }
  
  if (row_average) {
    penalty <- mean(penalty)
  }
  return(penalty)
}