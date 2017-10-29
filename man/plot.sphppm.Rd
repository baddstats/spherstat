\name{plot.sphppm}
\alias{plot.sphppm}
\title{
  Plot Fitted Point Process Intensity
}
\description{
  Plot the fitted intensity of a point process model
  as an image wrapped on the sphere.
}
\usage{
\method{plot}{sphppm}(x,
             eye = place("nedlands"), top = place("northpole"),
             w,
             \dots,
             eps = NULL, dimyx = NULL, main = "",
             action = c("image", "contour", "imagecontour"),
             col.image = NULL, col.lines = NULL)
}
\arguments{
  \item{x}{
    Fitted point process model on the sphere (object of class \code{"sphppm"}).
  }
  \item{eye}{
    Viewpoint. A vector of length 3 (or a pair of
    \code{list(lon,lat)} coordinates, in degrees)
    determining a position in space.
  }
  \item{top}{
    Vector of length 3 (or a pair of \code{list(lon,lat)} coordinates,
    in degrees) determining a position in space. The plot will be
    rotated so that this position appears to be directly above the
    centre of the earth. 
  }
  \item{w}{
    Spherical window (object of class \code{"sphwin"}) in which the
    fitted intensity should be plotted. Default is the window of the
    original data to which the model was fitted.
  }
  \item{\dots}{
    Additional arguments controlling the plot, passed
    to \code{\link[spatstat]{plot.im}} or
    \code{\link[spatstat]{contour.im}}.
  }
  \item{eps,dimyx}{
    Arguments passed to \code{\link[spatstat]{as.mask}}
    controlling the resolution of the final plotted image.
  }
  \item{main}{
    Main title. A character string, character vector, or expression.
  }
  \item{action}{
    Character string specifying how the fitted intensity should be
    plotted. See Details.
  }
  \item{col.image}{
    Optional. A vector of colour values (a character vector which can be
    interpreted as colours) or an object of class \code{"colourmap"}
    determining the colour map for the image plot, if
    \code{action="image"} or \code{action="imagecontour"}.
    Overrides the argument \code{col} if it is present.
  }
  \item{col.lines}{
    Optional. A single colour (a character string or integer
    which can be interpreted as a colour)
    to be used for the contour lines, if
    \code{action="contour"} or \code{action="imagecontour"}.
    Overrides the argument \code{col} if it is present.    
  }
}
\details{
  This function computes the intensity of the fitted model \code{x}
  at essentially every location on the sphere, and produces a
  colour image plot and/or contour plot, showing the sphere
  viewed from a particular spatial position.

  Options are \code{action="image"} to produce a colour image,
  \code{action="contour"} for a contour plot, and
  \code{action="imagecontour"} for both a colour image with contour
  lines superimposed.
  
  The algorithm runs as follows: First
  the intensity of the fitted model \code{x} is computed at
  a fine grid of locations on the sphere. Those grid locations which are
  visible from the viewing location \code{eye} are projected onto
  a two-dimensional plane. The values at these projected locations are
  interpolated onto a two-dimensional rectangular pixel grid by
  kernel smoothing, using \code{\link[spatstat]{Smooth.ppp}}, with a very small
  bandwidth \code{sigma}. The resulting two-dimensional pixel image is
  then plotted by \code{\link[spatstat]{plot.im}} and/or a contour map is drawn
  by \code{\link[spatstat]{contour.im}}.

  To control colours in the plot, the argument \code{col} can be used if
  this is unambiguous (i.e. if \code{action="image"} or
  \code{action="contour"}). Otherwise the arguments \code{col.image}
  and \code{col.lines} determine the colours used in the
  image and contour plots, respectively.

  Note that the plot shows the half of the sphere that is visible from the eye position (i.e. all locations \eqn{u} for which \eqn{d(u,eye) \leq \pi\rho/2}{d(u,eye) <= pi*rho/2}, where \eqn{\rho}{rho} is the radius of the sphere).
}
\value{
  Null.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{sphppm}},
  \code{\link{predict.sphppm}}
}
\examples{
  XX <- runif.sphwin(100)
  fit <- sphppm(XX ~ x1+ x2 + x3)
  plot(fit)
  plot(fit, dimyx=256)
  globedrawlong(seq(0, 350, by=10))
  globedrawlat(seq(-80,80,by=10))
  plot(fit, dimyx=256, col=heat.colors(128))
  plot(fit, dimyx=256, action="contour", col="green")
  plot(fit, dimyx=256, action="imagecontour",
                col.image=topo.colors(128), col.lines="white")
}
\keyword{spatial}
\keyword{hplot}
