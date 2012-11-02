#' @title Improved image() function
#' @aliases make_legend
#' @description 
#' Improved version of the \code{\link[graphics]{image}} function in 
#' the \code{graphics} package.
#' In particular, it displays matrices the way they are shown in the R console,
#' not transposed/rearranged/... For example, a covariance matrix has the
#' diagonal in from top-left to bottom-right as it should be, and not from 
#' bottom-left to top-right. 
#' 
#' The function \code{\link{make_legend}} also provides a better color scale 
#' legend handling.
#' 
#' Optionally \code{image2} displays a color histogram below the image,
#' which can be used to refine the display of a matrix by trimming outliers (as 
#' they can often distort the color representation).
#'
#' @param z a matrix containing the values to be plotted (NAs are allowed).
#' @param col colors: either a string decribing a pallette from the 
#' \code{RColorBrewer} package (see also \url{http://colorbrewer2.org/}), or a list of colors 
#' (see \code{\link[graphics]{image}} for suggestions).
#' @param axes a logical value indicating whether both axes should be 
#' drawn on the plot. 
#' @param xlab a label for the x axis 
#' @param ylab a label for the y axis
#' @param legend logical; if \code{TRUE} a color legend for will 
#' be plotted
#' @param zlim minimum and maximum z values for which colors should be 
#' plotted, defaulting to the range of the finite values of \code{z}.
#' @param zlim.label character string (default: \code{"color scale"}) to write 
#' next to the color legend
#' @param density logical; if \code{TRUE} a color histogram 
#' (\code{\link[stats]{density}}) will be plotted. Default: \code{FALSE}.
#' @param max.height height of the density plot (typically not modified by user)
#' @param ... optional arguments passed to \code{\link[graphics]{image}}
#' @keywords hplot aplot
#' @export
#' @seealso \code{\link[graphics]{image}}, \code{\link[fields]{image.plot}}
#' @examples
#' 
#' \dontrun{
#' # Correlation matrix has diagonal correctly
#' data(iris)
#' image(cor(as.matrix(iris[,-1])))
#' image2(cor(as.matrix(iris[,-1])), col = "RdBu")
#' }
#' # Color histogram
#' nn <- 100
#' set.seed(nn)
#' AA <- matrix(sample(c(rnorm(nn^2/2, -1, .1), rexp(nn^2/2, 1))), ncol = nn)
#' 
#' image2(AA, col = "Spectral", density = TRUE)
#' 
#' image2(AA, col = "Spectral", density = TRUE, zlim = c(min(AA), 3))
#' 
#' 

image2 <- function(z, col = NULL, 
                   axes = FALSE, legend = TRUE, xlab="", ylab = "", 
                   zlim = NULL, density = FALSE, max.height = NULL, 
                   zlim.label = "color scale", ...){
  
   zz <- z  
   zz[is.na(zz)] = min(zz, na.rm=TRUE)
   min <- min(zz)
   max <- max(zz)
   yLabels <- rownames(zz)
   xLabels <- colnames(zz)
   title <-c()
   if (!is.null(zlim)){
     min = zlim[1]
     max = zlim[2]
   }
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- pretty(1:ncol(zz))
  }
  if( is.null(yLabels) ){
    yLabels <- pretty(1:nrow(zz))
  }
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  if (is.null(col)) {
    col <- colorRampPalette(brewer.pal(9, name="YlOrRd"))(100)
  } 
  if (length(col) == 1 & is.character(col)){
    col <- colorRampPalette(brewer.pal(9, name=col))(100)
  }
  
  color_levels <- seq(min, max, length=length(col))

  # Reverse Y axis
  reverse <- nrow(zz) : 1
  yLabels <- yLabels[reverse]
  zz <- zz[reverse,]
  
  # Data Map
  if (legend) {
    if (density) {
      layout(matrix(data=c(1,2), nrow=2), heights=c(9,2))
    } else {
      layout(matrix(data=c(1,2), ncol=2), widths=c(9,2))
    }
    par(mar = c(1,1,2,1), cex.lab = 2, cex.axis = 2, lwd = 2)
  }
  
  # par(mar = rep(0.5, 4), cex.lab = 2, cex.axis = 2, lwd = 2)
  image(1:ncol(zz), 1:nrow(zz), t(zz), col=col, axes=FALSE, zlim=c(min,max), xlab = xlab, ylab = ylab, ...) 
  box()
  if( !is.null(title) ){
   title(main=title)
  }
  if (axes) {
    axis(1, at=xLabels, labels=xLabels, cex.axis=2)
    axis(2, at=yLabels, labels=yLabels, las=1, cex.axis=2)
  }
  
  
  if (legend) {
    # Color Scale
    if (density) {
      temp.pdf = density(zz)
      if (is.null(max.height)) {
        max.height = max(temp.pdf$y)*1.05
      }
      par(cex.lab = 2, cex.axis = 2, lwd = 2, mar = c(3,2.5,0,3))
      make_legend(data=zz, col = col, side = 1, zlim = c(min, max), 
                    col.ticks = pretty(temp.pdf$x, n = 6), cex.axis = 1.25, max.height = max.height )
      #legend("topright", col = "black", "frequency", lty = 1, lwd = 2, bty = "n", cex = 1.25)
      mtext(zlim.label, side = 1, line = 2, cex = 1.25)
      lines(temp.pdf, lwd = 2)
      axis(4, at = c(0, max.height/2, max.height), labels = c("0", paste(round(c(0.5, 1)* max.height, 1))) , cex.axis = 1.25)
      mtext("density", side = 2, line = 0.5, cex = 1.25)
      
      #mtext("bits", side = 1, line = 2, cex = 1.25)
      } else {
        par(cex.lab = 2, cex.axis = 2, lwd = 2, mar = c(1,1,2,4), las = 2)
        make_legend(data = zz, col = col, side = 4, zlim = zlim)
      }
  }
}

#' @rdname image2
#' @param data data for which the legend should be plotted
#' @param side on which side of the plot (1=bottom, 2=left, 3=top, 4=right)
#' @param col.ticks color tick marks
#' @param cex.axis The magnification to be used for axis annotation relative to 
#' the current setting of \code{cex}.
#' @param col.label same as \code{zlim.label}
#' @keywords aplot color

make_legend <- function(data = NULL, col = NULL, side = 1, 
                        zlim = NULL, col.ticks = NULL, cex.axis = 2, 
                        max.height = 1, col.label = "") {
  
  min <- min(data)
  max <- max(data)
  if (is.null(zlim)) {
    zlim <- c(min, max)
  }
  if (length(col) == 1 & is.character(col)) {
    col <- colorRampPalette(brewer.pal(9, name = col))(100)
  }
  color_levels <- seq(zlim[1], zlim[2], length = length(col))

  if (is.null(col.ticks)) {
    col.ticks <- pretty(color_levels)
    # col.ticks = c(floor(min*10)/10, pretty(color_levels), ceiling(max*10)/10)
    # col.ticks = unique(col.ticks)
  }
  
  if (any(side == c(1, 3))) {
    image(color_levels, c(0, max.height), 
          matrix(data = rep(color_levels, time = 2), nrow = length(color_levels), ncol = 2), 
          col = col, xlab = "", 
          ylab = "", axes = FALSE, zlim = zlim, ylim = c(0, max.height))
  } else {
    image(c(0, max.height), color_levels, 
          matrix(data = rep(color_levels, each = 2), ncol = length(color_levels), nrow = 2), 
          col = col, xlab = "", 
          ylab = "", axes = FALSE, zlim = zlim, xlim = c(0, max.height))
  }
  axis(side = side, at = col.ticks, labels = paste(col.ticks), cex.axis = cex.axis, 
       col = "black")
  mtext(col.label, side = side, line = 3, cex = 0.75 * cex.axis)
  # axis(side = side+1)
  box()
} 