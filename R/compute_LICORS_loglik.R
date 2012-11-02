#' @title Log-likelihood of LICORS model
#'
#' @description 
#' Computes the \emph{average} log-likelihood \eqn{\frac{1}{N} \ell(\mathbf{W}; \mathcal{D})} as a function of the 
#' weight matrix \eqn{\mathbf{W}} and the predictive state distributions 
#' \eqn{P(X = x \mid S = s_j) \approx P(X = x \mid \mathbf{W}_j)} for all 
#' \eqn{j = 1, \ldots, K}. See References.
#'
#' @param weight_matrix \eqn{N \times K} weight matrix
#' @param FLC_pdfs an \eqn{N \times K} matrix containing the estimates of all \eqn{K} 
#' FLC densities evaluated at all \eqn{N} sample FLCs.
#' @param lambda regularization parameter. Default: \code{lambda=0} 
#' (\code{penalty} and \code{q} will be ignored in this case).
#' @param penalty type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param base logarithm base for the \code{"entropy"} penalty.
#' Default: \code{base = 2}.  Any other real number is allowed; 
#' if \code{base = "nstates"} then it will internally assign it 
#' \code{base = ncol(weight_matrix)}.
#' @param q exponent for \eqn{L_q} norm.
#' @keywords manip nonparametric
#' @export
#' 


compute_LICORS_loglik <- function(weight_matrix = NULL, FLC_pdfs = NULL, 
                                  lambda = 0, penalty = "entropy", q = 2, base = exp(1)) {
  qq <- q
  if (is.null(weight_matrix) | is.null(FLC_pdfs) ){
    stop("You must provide the weight matrix and the FLC distributions.")
  }
  
  loglik <- mean(log(rowSums(FLC_pdfs * weight_matrix)))
  if (lambda > 0) {
    NN <- nrow(FLC_pdfs)
    kk <- ncol(FLC_pdfs)
    loglik <- loglik - lambda * compute_mixture_penalty(weight_matrix, 
                                                        type = penalty, 
                                                        q = qq, base = base) #- 1/2 * kk * log(NN)/NN
  }
  return(loglik)
} 