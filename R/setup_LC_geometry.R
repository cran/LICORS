#' @title Setup light cone geometry
#' @description 
#' \code{setup_LC_geometry} sets up the light cone geometry for LICORS.
#'
#' @param horizon a list with \code{PLC} and \code{FLC} horizon
#' @param speed speed of propagation
#' @param space_dim dimension of the spatial grid. Eg. \code{2} if the data is 
#' a video ( = image sequences).
#' @param shape shape of light cone: \code{'cone'}, \code{'tube'}, 
#' or \code{'revcone'}.
#' @keywords manip
#' @return 
#' A list of class \code{"LC"}.
#' 
#' @export
#' @seealso \code{\link{LC-utils}}, \code{\link{compute_LC_coords}}
#' @examples
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed=1, 
#'                        space_dim = 1, shape = "cone")
#' aa
#' plot(aa)
#' summary(aa)
#' 
setup_LC_geometry <- function(horizon = list(PLC = 1, FLC = 0), 
                              speed = 1, space_dim = 1, shape = "cone") {
  
  out <- list()
    
  if (is.null(horizon$FLC)) {
    horizon$FLC <- 0
  }
  
  out$horizon <- horizon
  out$speed <- speed
  out$space_dim <- 1
  out$shape <- "cone"
  
  out$coords <- list()
  out$coords$PLC <- compute_LC_coords(horizon = horizon$PLC, speed = speed, 
                                      space_dim = space_dim, shape = shape, 
                                      type = "PLC")
  out$coords$FLC <- compute_LC_coords(horizon = horizon$FLC, speed = speed, 
                                      space_dim = space_dim, shape = shape, 
                                      type = "FLC")
  out$n_p <- nrow(out$coords$PLC)
  out$n_f <- nrow(out$coords$FLC)
  class(out) <- "LC"
  return(out)
} 
