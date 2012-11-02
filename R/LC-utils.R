#' @title Utilities for LC class
#' @name LC-utils
#' @aliases plot.LC summary.LC LC_coords2control_settings
#' @description 
#' 
#' The \code{"LC"} class is the core property of the LICORS model as it specifies
#' the spatio-temporal neighborhood of the past and future light cone. 
#' The function \code{\link{setup_LC_geometry}} generates an instance of 
#' the \code{"LC"} class.
#' 
#' \code{plot.LC} plots LCs of \eqn{(1+1)D} and \eqn{(2+1)D} systems with 
#' respect to the origin \eqn{(\mathbf{r}, t) = (\boldsymbol 0, 0)}. 
#' This is especially useful for a quick check if the LC geometry 
#' specified by \code{\link{setup_LC_geometry}} is really the intended one.
#'
#' \code{summary.LC} prints a summary of the LC geometry.  
#' Returns (invisible) the summary matrix.
#'
#' \code{LC_coords2control_setting} computes auxiliary measures given the LC 
#' geometry such as horizon, spatial/temporal extension, etc. This function 
#' should not be called by the user directly; only by \code{\link{get_spacetime_iterator}}.
#' 
NULL

#' @rdname LC-utils
#' @keywords hplot
#' @method plot LC
#' @param x an object of class \code{"LC"} (see \code{\link{setup_LC_geometry}})
#' @param cex.axis The magnification to be used for axis annotation relative to 
#' the current setting of \code{cex}.
#' @param cex.lab The magnification to be used for x and y labels 
#' relative to the current setting of \code{cex}.
#' @param ... optional arguments passed to \code{plot}.
#' @export
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC= 1), speed=1, space_dim = 1, shape = "cone")
#' plot(aa)
#' bb = setup_LC_geometry(horizon = list(PLC = 2, FLC= 1), speed=1, space_dim = 1, shape = "revcone")
#' plot(bb)
#' 
plot.LC <- function(x, cex.axis = 2, cex.lab = 2, ...) {
  object <- x
  
  LC_coords <- rbind(object$coords$PLC, object$coords$FLC)
  if (object$space_dim == 0) {
    
  } else if (object$space_dim == 1) {
    par(mar = c(4, 5, 1, 1))
    plot(LC_coords, type = "n", cex.lab = cex.lab, xlab = "Time", ylab = "Space", 
         axes = FALSE, ...)
    box()
    abline(v = 0)
    
    axis(1, at = unique(LC_coords[, 1]), cex.axis = cex.axis)
    axis(2, at = unique(LC_coords[, 2]), cex.axis = cex.axis)
    grid()
    points(object$coords$FLC, pch = "+", cex = 3, col = "blue")
    points(object$coords$PLC, pch = "--", cex = 3, col = "red")
    text(-0.75, max(LC_coords[, 2]) - 0.5, "PLC", col = "red", cex = 2)
    text(-0.25, min(LC_coords[, 2]) + 0.5, "FLC", col = "blue", cex = 2)
  } else if (object$space_dim == 2) {
    print("Not yet implemented.")
  } else {
    plot.new()
    legend("top", paste("Can't plot (3+1)D systems \n not possible.", cex = 2))
  }
  
  # if (object$space_dim < 3){ title('PLC (red) and FLC (blue) shape') }
} 

#' @rdname LC-utils
#' @keywords models print
#' @param object an object of class \code{"LC"}
#' @param verbose logical; if \code{TRUE} LC information is printed in the 
#' console
#' @export
#' @method summary LC
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed=1, space_dim = 1, shape = "cone")
#' summary(aa)
#' 
summary.LC <- function(object, verbose = TRUE, ...) {
  
  temp_mat <- rbind(object$horizon)
  temp_mat <- rbind(temp_mat, object$speed)
  temp_mat <- as.data.frame(temp_mat)
  temp_mat <- rbind(temp_mat, object$shape)
  temp_mat <- rbind(temp_mat, c(object$n_p, object$n_f))
  
  colnames(temp_mat) <- c("PLC", "FLC")
  rownames(temp_mat) <- c("horizon", "speed", "shape", "dimensionality")
  
  if (verbose) {
    cat(rep("*", 20))
    cat("\n \n")
    cat(paste("The field extends over a ", object$space_dim, "-dimensional space. \n", 
              sep = ""))
    cat("Light cones have therefore the following characteristics: \n \n")
    print(temp_mat)
    cat("\n")
    cat(rep("*", 20))
  }
  
  out <- list()
  out$table <- temp_mat
  class(out) <- "summary.LC"
  invisible(out)
}

#' @rdname LC-utils
#' @param LC_coords template of a light cone (with respect to origin)
#' @keywords method
#' @export
#' @seealso \code{\link{compute_LC_coords}}
#' @examples
#' 
#' aa = setup_LC_geometry(horizon = list(PLC = 2, FLC = 0), speed=1, space_dim = 1, shape = "cone")
#' LC_coords2control_settings(aa$coords)
#' 

LC_coords2control_settings <- function(LC_coords) {
  
  out <- list()
  
  out$horizon <- list()
  out$horizon$PLC <- -min(LC_coords$PLC[, 1])
  out$horizon$FLC <- -min(LC_coords$FLC[, 1])
  
  LCs_coords <- rbind(LC_coords$PLC, LC_coords$FLC)
  out$space_dim <- ncol(LCs_coords) - 1
  if (out$space_dim == 1) {
    out$space_cutoff <- max(LCs_coords[, -1])
  } else {
    out$space_cutoff <- apply(LCs_coords[, -1], 2, max)
  }
  out$n_p <- nrow(LC_coords$PLC)
  out$n_f <- nrow(LC_coords$FLC)
  
  return(out)
  
}

