\name{image_scale}
\alias{image_scale}
\alias{image.scale}
\encoding{UTF-8}
\title{
Helper function to draw scale for image
}
\description{
Plots a scale in the right margin of a plot, typically an image plot.
}
\usage{
image_scale(z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", "ranges"))
}
\arguments{
  \item{z}{
the values plotted in the image, the \code{z} argument to \code{image}
}
  \item{col}{
a vector of colours, use the same as in the image
}
  \item{x}{
the x coordinate of the left edge of the scale bar; or a list with components x and y, each of length 2, giving the x and y coordinates of the edges of the scale bar.
}
  \item{y}{
the y coordinate of the bottom edge of the scale bar, ignored if \code{x} is a list
}
  \item{size}{
the size of the boxes making up the scale bar, a length 2 vector with the width and height, or a scalar, in which case width = height; ignored if \code{x} is a list
}
  \item{digits}{
the number of decimal places to display
}
  \item{labels}{
if "breaks", the dividing lines in the scale bar are labelled, if "ranges", each box is labelled with its range; may be abbreviated
}
}

\value{
None, used for its plotting side effect.
}
\references{
Royle, J.A., Chandler, R.B., Sollmann, R., & Gardner, B. (2014) \emph{Spatial capture-recapture} Elsevier.
}
\author{
Borrowed from the \pkg{SCRbook} package.
}
\note{
This function appears in the text as \code{image.scale}; renamed to avoid confusion with generic \code{image} functions.
}

\examples{
# uses the built-in volcano data set
require(grDevices) # for colours
require(graphics)
par(mar = c(3,3,3,6))  # make the right margin wide enough
image(t(volcano)[ncol(volcano):1,], col=terrain.colors(12)) 
image_scale(volcano, col=terrain.colors(12))
# Try placing the scale bar on the left
par(mar = c(3,8,3,1))  # make the left margin wide enough
image(t(volcano)[ncol(volcano):1,], col=terrain.colors(12)) 
image_scale(volcano, col=terrain.colors(12), x= -0.28, digits=0)
# Trial and error needed to get the x value right.
}
