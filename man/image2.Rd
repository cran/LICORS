\name{image2}
\alias{image2}
\alias{make_legend}
\title{Improved image() function}
\usage{
image2(x = NULL, y = NULL, z = NULL, col = NULL, axes = FALSE, legend = TRUE, 
    xlab = "", ylab = "", zlim = NULL, density = FALSE, max.height = NULL, 
    zlim.label = "color scale", ...)

make_legend(data = NULL, col = NULL, side = 1, zlim = NULL, col.ticks = NULL, 
    cex.axis = 2, max.height = 1, col.label = "")
}
\arguments{
  \item{x,y}{locations of grid lines at which the values in
  z are measured. These must be finite, non-missing and in
  (strictly) ascending order. By default, equally spaced
  values from 0 to 1 are used. If x is a list, its
  components \code{x$x} and \code{x$y} are used for
  \code{x} and \code{y}, respectively. If the list has
  component \code{z} this is used for \code{z}.}

  \item{z}{a matrix containing the values to be plotted
  (NAs are allowed). Note that \code{x} can be used instead
  of \code{z} for convenience.}

  \item{col}{colors: either a string decribing a pallette
  from the \code{RColorBrewer} package (see also
  \url{http://colorbrewer2.org/}), or a list of colors (see
  \code{\link[graphics]{image}} for suggestions).}

  \item{axes}{a logical value indicating whether both axes
  should be drawn on the plot.}

  \item{xlab}{a label for the x axis}

  \item{ylab}{a label for the y axis}

  \item{legend}{logical; if \code{TRUE} a color legend for
  will be plotted}

  \item{zlim}{minimum and maximum z values for which colors
  should be plotted, defaulting to the range of the finite
  values of \code{z}.}

  \item{zlim.label}{character string (default: \code{"color
  scale"}) to write next to the color legend}

  \item{density}{logical; if \code{TRUE} a color histogram
  (\code{\link[stats]{density}}) will be plotted. Default:
  \code{FALSE}.}

  \item{max.height}{height of the density plot (typically
  not modified by user)}

  \item{...}{optional arguments passed to
  \code{\link[graphics]{image}}}

  \item{data}{data for which the legend should be plotted}

  \item{side}{on which side of the plot (1=bottom, 2=left,
  3=top, 4=right)}

  \item{col.ticks}{color tick marks}

  \item{cex.axis}{The magnification to be used for axis
  annotation relative to the current setting of
  \code{cex}.}

  \item{col.label}{same as \code{zlim.label}}
}
\description{
  Improved version of the \code{\link[graphics]{image}}
  function in the \code{graphics} package. In particular,
  it displays matrices the way they are shown in the R
  console, not transposed/rearranged/... For example, a
  covariance matrix has the diagonal in from top-left to
  bottom-right as it should be, and not from bottom-left to
  top-right.

  The function \code{\link{make_legend}} also provides a
  better color scale legend handling.

  Optionally \code{image2} displays a color histogram below
  the image, which can be used to refine the display of a
  matrix by trimming outliers (as they can often distort
  the color representation).
}
\examples{
\dontrun{
# Correlation matrix
data(iris)  # make sure its from 'datasets' package, not from 'locfit'
image(cor(as.matrix(iris[, names(iris) != "Species"])))

# Correlation matrix has diagonal from top left to bottom right
par(mar = c(1, 3, 1, 2))
image2(cor(as.matrix(iris[, names(iris) != "Species"])), col = "RdBu", axes = FALSE)
}
# Color histogram
nn <- 10
set.seed(nn)
AA <- matrix(sample(c(rnorm(nn^2, -1, 0.1), rexp(nn^2/2, 0.5))), ncol = nn)

image2(AA, col = "Spectral")
image2(y = 1:15 + 2, x = 1:10, AA, col = "Spectral", axes = TRUE)
image2(y = 1:15 + 2, x = 1:10, AA, col = "Spectral", density = TRUE, axes = TRUE)

image2(AA, col = "Spectral", density = TRUE, zlim = c(min(AA), 3))
}
\seealso{
  \code{\link[graphics]{image}},
  \code{\link[fields]{image.plot}}
}
\keyword{aplot}
\keyword{color}
\keyword{hplot}
