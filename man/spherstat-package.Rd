\name{spherstat-package}
\alias{spherstat-package}
\alias{spherstat}
\docType{package}
\title{
  Analysis of Point Patterns on a Sphere
}
\description{
  A package for analysing spatial point pattern data on the sphere,
  or observed in a region of the sphere. 
}
\details{
  This package allows a user to analyse point pattern data on a sphere,
  or a subset of a sphere.  Subsets of the sphere that are currently
  supported are the spherical cap, band, wedge, and polygon.  Currently
  this package supports:
  \itemize{
    \item calculations in spherical geometry and trigonometry;
    \item calculation of pairwise and nearest neighbour distances
    between points;
    \item calculation of distance from point to boundary of window;
    \item estimation of summary functions including
    the spherical analogues of Ripley's \eqn{K} function,
    the nearest neighbour distance distribution function \eqn{G},
    the empty space function \eqn{F} and the \eqn{J} function;
    \item
    simulation of Poisson, Matern inhibition, Strauss, Hard-core and Neyman-Scott point processes;
    \item
    model-fitting for Poisson point processes and
    Neyman-Scott processes.
  }
  The preferred coordinate system for data given to functions in this
  package is the spherical coordinate system.  In other words, the
  coordinate of a point should be given as \eqn{(\theta, \phi)}{(theta, phi)}, where
  \eqn{\theta}{theta} is the colatitude, between \code{0} and \eqn{\pi}{pi},
  and \eqn{\phi}{phi} is longitude and is between \code{0} and \eqn{2\pi}{2 * pi}.
  The north pole has colatitude \code{0} and arbitrary
  longitude (although for most cases, use \eqn{\phi=0}{phi=0}), and longitudes are
  measured in anticlockwise direction.  Some functions allow data to be
  given using Cartesian coordinates \code{(x,y,z)}.
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
Adrian Baddeley
}
\note{
  Many functions in these package are the analogues for point processes on
  the sphere of functions available in the package \pkg{spatstat} for
  point processes in \eqn{R^2}{R^2} or \eqn{R^3}{R^3}.  Elements
  of the code and help pages in this package have been taken
  with permission from
  corresponding material in \pkg{spatstat}.
}
\keyword{spatial}
\examples{
  sph <- sphwin(type="sphere")
  sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
  sph.Gest <- Gsphere(X=sph.pp, win=sph)
}
