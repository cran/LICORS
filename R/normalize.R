#' @title Normalize a matrix/vector to sum to one (probability simplex)
#'
#' @description 
#' \code{normalize} projects a vector or matrix onto the probability simplex.
#' 
#' If all entries (per row or column) get thresholded to \eqn{0} (since they are 
#' all negative to start with), then it 
#' sets the position of the maximum of \code{x} to \eqn{1} and leaves all other
#' entries at \eqn{0}.
#'
#' @param x a numeric matrix(like object).
#' @param byrow logical; if \code{TRUE} rows are normalized; otherwise columns.
#' @param tol a tolerance level to set values \eqn{< tol} to \eqn{0} (after
#' an initial normalization). Default: \code{tol=1e-6}
#' @keywords manip array
#' @export
#' @seealso \code{\link{threshold}}
#' @return
#' If \code{x} is a vector it returns the thresholded vector 
#' (see \code{\link{threshold}}) and normalized by its sum.
#' If \code{x} is a matrix it works by column of by row 
#' (argument \code{byrow}).
#' @examples
#' print(normalize(c(1,4,2,2,10)))
#' print(normalize(c(-1,-2, -1)))
#' AA = matrix(rnorm(12), ncol = 3)
#' print(normalize(AA, byrow = TRUE))
#' print(normalize(AA, byrow = FALSE))

normalize <- function(x, byrow = TRUE, tol = 1e-6) {
  
  object <- x
  
  if (byrow == FALSE){
    return( t(normalize( t(object), tol = tol, byrow = TRUE) ) )
  } else {
    if (is.matrix(object)){
      matrix.new <- threshold(object, min = 0)
      max.pos <- integer(0)
      if (any(object < 0)){
        all.zeros <- which(apply(matrix.new, 1, function(u) all(u == 0)))
        if (length(all.zeros) > 1){
          max.pos <- apply(object[all.zeros,], 1, which.max)
        } else if (length(all.zeros) == 1) {
          max.pos = which.max(object[all.zeros,])
        }
        matrix.new[cbind(all.zeros, max.pos)] = 1
      }
      
      #  if (any(all.zeros)) {
      #      matrix.new[all.zeros, ] <- runif(sum(all.zeros) * ncol(matrix.new), 
      #          0, 1)
      #  }
      matrix.new <- sweep(matrix.new, 1, rowSums(matrix.new), "/")
      if (tol > 0) {
        matrix.new[matrix.new < tol] <- 0
        matrix.new <- sweep(matrix.new, 1, rowSums(matrix.new), "/")
      } 
      invisible(matrix.new)
    } else if (is.vector(object)){
      max.pos = which.max(object)
      max.val = object[max.pos]
      object = threshold(object, min = 0)
      if (all(object == 0)){
        object[max.pos] = 1
      } else {
        object = object/sum(object)
      }
      invisible(object)
    }
  }
} 
