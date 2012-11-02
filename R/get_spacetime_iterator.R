#' @title 
#' Get an iterator over the space-time coordinates of the field
#' @description 
#' This function returns iterators over space-time coordinates of the field.
#' Both for the whole field as well as the truncated field (without the margin)
#' @param dim dimension of the original field (first dimension is time; rest is 
#' space)
#' @param LC_coords template of the LC coords
#' @keywords iteration
#' @export
#' @seealso \code{\link{compute_LC_coords}}, \code{\link{setup_LC_geometry}}
#' @import itertools
#' @examples
#' AA = matrix(rnorm(200), ncol = 10)
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 3, FLC = 0), shape = "cone")
#' bb = get_spacetime_iterator(dim(AA), LC_geom$coords)
#' 

get_spacetime_iterator <- function(dim = NULL, LC_coords = NULL) {
  if (is.null(dim)) {
    stop("You must provide the dimensions of the original field.")
  }
  
  controls <- LC_coords2control_settings(LC_coords)
  horizon <- controls$horizon
  
  space_cardinality <- length(dim) - 1
  TT <- dim[1]
  
  space_dim <- list()
  for (ss in 1:space_cardinality) {
    space_dim[[ss]] <- dim[-1][ss]
  }
  
  truncated_dim <- compute_margin_coords(dim, LC_coords)
  
  if (space_cardinality == 1) {
    iterator_spacetime <- ihasNext(product(tt = 1:TT, x1 = 1:space_dim[[1]]))
    iterator_spacetime.truncated <- ihasNext(product(tt = horizon$PLC + 1:(TT - 
      sum(unlist(horizon))), x1 = truncated_dim$space[[1]][1]:truncated_dim$space[[1]][2]))
  } else if (space_cardinality == 2) {
    iterator_spacetime <- ihasNext(product(tt = 1:TT, x1 = 1:space_dim[[1]], 
                                           x2 = 1:space_dim[[2]]))
    iterator_spacetime.truncated <- ihasNext(product(tt = horizon$PLC + 1:(TT - 
      sum(unlist(horizon))), x1 = truncated_dim$space[[1]][1]:truncated_dim$space[[1]][2], 
                                                     x1 = truncated_dim$space[[2]][1]:truncated_dim$space[[2]][2]))
  } else {
    stop("Space dimension > 2 is not implemented.")
  }
  out <- list()
  out$all <- iterator_spacetime
  out$dim_all <- dim
  out$length_all <- prod(out$dim_all)
  out$truncated <- iterator_spacetime.truncated
  out$dim_truncated <- truncated_dim$dim
  out$length_truncated <- prod(out$dim_truncated)
  
  names(out$dim_all) = names(out$dim_truncated) = c("time", paste("space", 1:length(space_dim), sep=""))
  
  return(out)
}