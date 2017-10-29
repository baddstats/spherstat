\name{plot.sphwin}
\alias{plot.sphwin}
\title{
  Plot a Spherical Window
}
\description{
  Plots a spherical window, representing a region on the sphere.
}
\usage{
 \method{plot}{sphwin}(x, \dots, eye, top, add = FALSE)
}
\arguments{
  \item{x}{
    Object of class \code{"sphwin"} representing a region on the sphere.
  }
  \item{\dots}{
    Additional graphics parameters.
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
    Logical value indicating whether to add the window to an
    existing plot. This will only work correctly if the same
    values of \code{eye} and \code{top} are used in both plots.
  }
}
\details{
  The plot shows the half of the sphere that is visible
  from the \code{eye} position (i.e. all locations \eqn{u} for which \eqn{d(u,eye) \leq \pi\rho/2}{d(u,eye) <= pi*rho/2}, where \eqn{\rho}{rho} is the radius of the sphere).
}
\value{
  Null.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{plot.sp2}}, 
  \code{\link{sphwin}}.
}
\examples{
  sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(pi/8,pi/6))
  plot(sphcap)
}
\keyword{spatial}
\keyword{hplot}

