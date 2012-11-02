#' @title Computes coordinates of PLC and FLC relative to origin
#'
#' @description 
#' Computes the space-time coordinates of PLC and FLC given control settings
#' relative to the origin \eqn{(\mathbf{r}, t) = (\boldsymbol 0, 0)}. 
#' 
#' Since these coordinates do not change for different space-time positions, 
#' they can be computed once before getting the LC configurations for the entire
#' field and then used in each call by array indexing in 
#' \code{\link{get_LC_config}}.
#'
#'
#' @param horizon integer; horizon for the PLC or FLC
#' @param speed speed of propagation
#' @param space_dim maximum value 
#' @param type \code{"PLC"} or \code{"FLC"}
#' @param shape shape of light cone: \code{'cone'}, \code{'tube'}, 
#' or \code{'revcone'}.
#' @keywords manip
#' @export
#' @seealso \code{\link{get_LC_config}} \code{\link{setup_LC_geometry}} 
#' \code{\link{summary.LC}} \code{\link{plot.LC}}
#' 
#' @examples
#' plot(compute_LC_coords(speed=1, horizon = 4), xlim = c(-4, 2), pch = "-", 
#'      cex=2, col = 2,xlab = "Time", ylab = "Space")
#' points(compute_LC_coords(speed=1, horizon = 2, type = "FLC"), pch = "+", 
#'        cex=2, col ="blue")
#' 
#' plot(compute_LC_coords(speed=1, horizon = 4, shape = "tube", type = "FLC"))
#' plot(compute_LC_coords(speed=1, horizon = 4, shape = "revcone", type = "PLC"))
#' 

compute_LC_coords = function(horizon=1, speed=1, space_dim=1, type = "PLC", 
                          shape = "cone"){
  ind = matrix(0, ncol = space_dim + 1, nrow = 1)
  if (horizon == 0){
    return(ind)
  }
  if (shape == "cone"){
    for (tt in 1:horizon){
      ind = rbind(ind, cbind(-tt, seq(-tt*speed, tt*speed, 1)))
    }
  } else if (shape == "tube") {
    for (tt in 1:horizon){
      ind = rbind(ind, cbind(-tt, seq(-speed, speed, 1)))
    } 
  } else if (shape == "revcone"){
    for (tt in 0:horizon){
      ind = rbind(ind, cbind(-tt, seq(-(horizon-tt)*speed, (horizon-tt)*speed, 1)))
    }
    if (type == "PLC"){
      ind[,1] = ind[,1]-1
    }
    
  }
  
  if (type == "FLC"){
    ind[,1] = -ind[,1]
    return(ind)
  } else {
    return(ind[ind[,1]<0 & ind[,1] >= -horizon,])
  }
}