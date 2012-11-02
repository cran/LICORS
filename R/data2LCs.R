#' @title Iterate over (N+1)D field and get all LC configurations
#' 
#' @description 
#' \code{data2LCs} gets all PLC or FLC configuration from a \eqn{(N+1)D} field 
#' given the LC template.  The shape and dimension of this LC template 
#' depends on coordinates passed on by \code{\link{setup_LC_geometry}}.
#' 
#' \subsection{User-defined LC template}{
#' 
#' Since \code{data2LCs} passes the \code{LC_coords} array to 
#' \code{\link{get_LC_config}} to iterate over the entire dataset, this 
#' functional programming approach allows user-defined light cone shapes 
#' (independent of the shapes implemented by \code{\link{setup_LC_geometry}}). 
#' 
#' Just replace the \code{$coords} from the \code{"LC"} class with a
#' user-specified LC template.
#' }
#' 
#' @param field spatio-temporal field; either a matrix or a 3-dimensional array 
#' with time \eqn{t} as the first dimension, and the spatial coordinates as 
#' subsequent dimensions. Make sure to check \code{\link{compute_LC_coords}} 
#' for correct formatting.
#' @param LC_coords coordinates for LC shape and dimension (usually the 
#' \code{$coords} value from the \code{"LC"} class; but also user-defined
#' coordinates are possible here).
#' @keywords manip
#' @import itertools
#' @export
#' @seealso \code{\link{compute_LC_coords}}, \code{\link{setup_LC_geometry}}
#' @examples
#' set.seed(1)
#' AA = matrix(rnorm(200), ncol = 10)
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 2, FLC = 0), shape ="cone")
#' bb = data2LCs(AA, LC_coords = LC_geom$coords)
#' image2(bb$PLC)
#' plot(density(bb$FLC))
#' 
data2LCs <- function(field = NULL, 
                     LC_coords = list(PLC = NULL, FLC = NULL)){

  controls = LC_coords2control_settings(LC_coords)
  n_p = controls$n_p
  n_f = controls$n_f
  space_dim = controls$space_dim
  
  iter_spacetime = get_spacetime_iterator(dim(t(field)), LC_coords)
  
  LC_array = matrix(NA, nrow = iter_spacetime$length_truncated,
                    ncol = 1 + 1 + space_dim + n_f + n_p )
  colnames(LC_array) = c("ID", "time", 
                         paste("x", 1:space_dim, sep=""), 
                         paste("FLC",1:n_f, sep=""), 
                         paste("PLC", 1:n_p, sep=""))

  ii = 0
  
  if (controls$n_f == 1) {
    while (hasNext(iter_spacetime$truncated)){
      temp_coord = unlist(nextElem(iter_spacetime$truncated))
      ii = ii + 1
      LC_array[ii,1] = ii
      LC_array[ii, 1 + 1:length(temp_coord)] = temp_coord
      LC_array[ii, 1 + space_dim+1 + 1:n_f] = c(field[temp_coord[2], temp_coord[1]])
      LC_array[ii, 1 + space_dim+1 + n_f + 1:n_p] = get_LC_config(temp_coord, field, LC_coords$PLC)
    }
  } else {
    while (hasNext(iter_spacetime$truncated)){
      temp_coord = unlist(nextElem(iter_spacetime$truncated))
      ii = ii + 1
      LC_array[ii, 1] = ii
      LC_array[ii, 1 + 1:length(temp_coord)] = temp_coord
      LC_array[ii, 1 + space_dim+1 + 1:n_f] = get_LC_config(temp_coord, field, LC_coords$FLC) #c(field[temp_coord[2], temp_coord[1]])
      LC_array[ii, 1 + space_dim+1 + n_f + 1:n_p] = get_LC_config(temp_coord, field, LC_coords$PLC)
    }
  }
  out = list()
  out$FLC = cbind(LC_array[, 1 + space_dim+1 + 1:n_f])
  out$PLC = LC_array[, 1 + space_dim+1 + n_f + 1:n_p]
  #out$coords = LC_array[, c(1, 1+ 1:length(temp_coord))]
  
  invisible( out )
}





