#' @title Converts label vector to 0/1 mixture weight matrix 
#' 
#' @description 
#' Converts unique cluster assignment stored in a length \eqn{N} label vector 
#' into a \eqn{N \times K} Boolean matrix of mixture weights.
#' 
#' @param state_vector a vector of length \eqn{N} with the state labels
#' @param nstates_total total number of states. If \code{NULL}, then the 
#' maximum of \code{state_vector} is chosen
#' @keywords manip array
#' @export
#' @seealso \code{\link{weight_matrix2states}}
#' @examples
#' ss = sample.int(5, 10, replace = TRUE)
#' WW = states2weight_matrix(ss)
#' 
#' image2(WW, col = "RdBu", xlab = "States", ylab = "Samples", axes = FALSE)
#'       

states2weight_matrix <- function(state_vector=NULL, nstates_total = NULL){
  if (is.null(nstates_total)){
    nstates_total = max(state_vector)
  }
  weight_matrix = matrix(0, nrow = length(state_vector), ncol = nstates_total)
  weight_matrix[cbind(1:length(state_vector), state_vector)] = 1
  
  invisible(weight_matrix)
}

