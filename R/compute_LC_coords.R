#' @title Computes coordinates of PLC and FLC relative to origin
#'
#' @description 
#' Computes the space-time coordinates of PLC and FLC given control settings
#' relative to the origin \eqn{(\mathbf{r}, t) = (\boldsymbol 0, 0)}. 
#' 
#' Since these coordinates do not change for different space-time positions, 
#' they can be computed once before getting the LC configurations for the entire
#' field and then used in each call by array maskexing in 
#' \code{\link{get_LC_config}}.
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

compute_LC_coords <- function(horizon = 1, speed = 1, 
                              space_dim = 1, 
                              type = c("PLC", "FLC"),
                              shape = c("cone", "tube", "revcone")){
  
  if (!(space_dim %in% c(1, 2))) {
    stop("Space dimension ", space_dim, " is not supported yet.")
  }
  
  type <- match.arg(type)
  shape <- match.arg(shape)
  
  mask <- matrix(0, ncol = space_dim + 1, nrow = 1)
  if (horizon == 0){
    return(mask)
  }
  
  if (space_dim == 1) {
    switch(shape,
           cone = {
             for (tt in 1:horizon){
               mask <- rbind(mask,
                             cbind(-tt, seq(-tt*speed, tt*speed, 1))
                             )
             }
           },
           tube = {
             for (tt in 1:horizon){
               mask <- rbind(mask, 
                             cbind(-tt, seq(-speed, speed, 1))
                             )
             }
           },
           revcone = {
             for (tt in 0:horizon){
               mask <- rbind(mask, 
                             cbind(-tt, 
                                  seq(-(horizon-tt)*speed, (horizon-tt)*speed, 1))
                            )
             }
             if (type == "PLC"){
               mask[, 1] <- mask[, 1] - 1
             }
           }
           )
  } else { # for space_dim = 2
    switch(shape,
           cone = {       
             for (tt in seq_len(horizon)) {
               space.mask <- expand.grid(seq(-tt*speed, tt*speed, 1), 
                                         seq(-tt*speed, tt*speed, 1))
               mask <- rbind(mask, cbind(-tt, as.matrix(space.mask)))
             }
           },
           tube = {
             space.mask <- expand.grid(seq(-speed, speed, 1), 
                                       seq(-speed, speed, 1))
             for (tt in seq_len(horizon)) {
               mask <- rbind(mask, cbind(-tt, as.matrix(space.mask)))
             }
           },
           revcone = {
             for (tt in 0:horizon) {
               space.mask <- expand.grid(seq(-(horizon-tt)*speed, (horizon-tt)*speed, 1), 
                                         seq(-(horizon-tt)*speed, (horizon-tt)*speed, 1))
               mask <- rbind(mask, cbind(-tt, as.matrix(space.mask)))
             }
             if (type == "PLC"){
               mask[, 1] <- mask[, 1] - 1
             }
           }
    )
  }
  colnames(mask) <- c("time", paste0("x", seq_len(space_dim)))
  if (type == "FLC"){
    mask[, 'time'] <- -mask[, 'time']
    return(mask)
  } else {
    return(mask[mask[, 'time'] < 0 & mask[, 'time'] >= -horizon, ])
  }
}