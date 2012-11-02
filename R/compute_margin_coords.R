#' @title Get LC configuration from a (N+1)D field
#'
#' @description 
#' \code{compute_margin_coords} computes the coordinates (boundary) of the 
#' margin of the field.
#' 
#' @param dim a vector with the dimensions of the field (time, space1, space2, ..., spaceN)
#' @param LC_coords template of the LC coords
#' @keywords manip
#' @export
#' @seealso \code{\link{compute_LC_coords}}
#' @examples
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 3, FLC = 0), shape = "cone")
#' 
#' data(contCA00)
#' 
#' aa = compute_margin_coords(dim(contCA00$observed), LC_geom$coords)
#' aa
compute_margin_coords <- function(dim = NULL, LC_coords = NULL){
  if (is.null(dim) | is.null(LC_coords)){
    stop("You must provide dimensions of the field, and control settings speed & horizon.")
  }
  
  controls = LC_coords2control_settings(LC_coords)
  horizon = c(controls$horizon$PLC, controls$horizon$FLC)
  if (length(horizon) == 1){
    horizon = c(horizon, 0)
  }
  
  space_cardinality = length(dim) - 1
  TT = dim[1]
  space_dim = list()
  for (ss in 1:space_cardinality){
    space_dim[[ss]] = dim[-1][ss]
  }
  
  out = list()
  out$time = c(horizon[1] + 1, TT - horizon[2])
  out$space = list()
  out$dim = c(diff(out$time)+1)
  for (ss in 1:space_cardinality){
    out$space[[ss]] = c(controls$space_cutoff+1, space_dim[[ss]] - controls$space_cutoff )
    out$dim = c(out$dim, diff(out$space[[ss]]) + 1)
  }
  
  return(out)
}


