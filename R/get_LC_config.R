#' @title Get configuration of a light cone (LC)
#'
#' @description 
#' \code{get_LC_config} obtains the PLC or FLC at a particular 
#' \eqn{(\mathbf{r}, t)} from a \eqn{(N+1)D} field based on the LC 
#' template from \code{\link{compute_LC_coords}} (or 
#' \code{\link{setup_LC_geometry}}).
#' 
#' @param coord space-time coordinate \eqn{(\mathbf{r}, t)}
#' @param field spatio-temporal field; either a matrix or a 3-dimensional array 
#' with time \eqn{t} as the first coord, and the spatial coords in order.  
#' Make sure to see also \code{\link{compute_LC_coords}} for correct 
#' formatting.
#' @param LC_coords template coords for the LC
#' @keywords method utilities
#' @export
#' @seealso \code{\link{compute_LC_coords}}
#' @examples
#' AA = matrix(rnorm(40), ncol = 5)
#' image2(AA)
#' LCind = compute_LC_coords(speed = 1, horizon = 2, shape = "cone")
#' AA
#' get_LC_config(c(2,5), AA, LCind)
#' 
get_LC_config <- function(coord, 
                          field = NULL, 
                          LC_coords = NULL) {
  
  return( field[sweep(LC_coords[,2:1], 2, rbind(coord[2:1]), "+")] )
}

