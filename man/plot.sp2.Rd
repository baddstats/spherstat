\name{plot.sp2}
\alias{plot.sp2}
\alias{plot.sp3}
\title{
  Plot a Point Pattern on the Sphere
}
\description{
  Plot a point pattern on the sphere, and its window of observation.
}
\usage{
 \method{plot}{sp2}(x, \dots, eye, top, add = FALSE)

 \method{plot}{sp3}(x, \dots, eye, top, add = FALSE)
}
\arguments{
  \item{x}{
    Point pattern on the sphere (object of class \code{"sp2"} or \code{"sp3"}).
  }
  \item{\dots}{
    Graphical parameters passed to lower-level plot functions.
  }
  \item{eye}{
    Spatial position of the viewer. Either a numeric vector
    of length 2 (interpreted as giving colatitude and longitude
    respectively in radians),
    or a numeric vector of length 3 (interpreted as the
    three-dimensional Cartesian coordinates) or a list with entries
    \code{lon} and \code{lat} giving longitude and latitude coordinates
    in degrees.
  }
  \item{top}{
    Location on the sphere that should be at the top of the plot.
    Same format as \code{eye}. Defaults to the North Pole.
  }
  \item{add}{
    Logical value indicating whether to add the points to an
    existing spherical plot.
    This will only work correctly if the same
    values of \code{eye} and \code{top} are used in both plots.
  }
}
\details{
  Only points that are visible from the \code{eye} position
  are plotted (i.e. all points \eqn{x} for which \eqn{d(x,eye) \leq \pi\rho/2}{d(x,eye) <= pi*rho/2}, where \eqn{\rho}{rho} is the radius of the sphere). That is, points on the other side of the sphere
  are not displayed (i.e. all points \eqn{y} for which \eqn{d(y,eye) > \pi\rho/2}{d(x,eye) > pi*rho/2}).
}
\value{
  Null.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{plot.sphwin}} for plotting the window only.

  \code{\link{sp2}}, \code{\link{sp3}}
}
\examples{
  sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(pi/8,pi/6))
  X <- runif.sphwin(n=150, win=sphcap, as.sp=TRUE, ndim="2")
  plot(X)

  plot(X$win, col="pink")
  plot(X, add=TRUE, pch=".", cex=4, col="blue")
}
\keyword{spatial}
\keyword{hplot}

