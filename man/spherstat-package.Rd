\name{spherstat-package}
\alias{spherstat-package}
\alias{spherstat}
\docType{package}
\title{
Analysis of spatial point pattern data on spheres and subsets of spheres
}
\description{
This package allows a user to analyse point pattern data on a sphere, or a subset of a sphere.  Subsets of the sphere that are currently supported are the spherical cap, band, wedge, and polygon.  Currently this package allows pairwise, boundary and nearest neighbour distances to be calculated, Ripley's K, the nearest neighbour (G) and empty space (F) functions to be estimated and various point processes and densities on the sphere simulated.  In addition, some spherical geometric and trigonometric functions are available.

The preferred coordinate system for data given to functions in this package is the spherical coordinate system.  In other words, the coordinate of a point should be given as (theta, phi), where theta is the colatitude is between 0 and pi, and phi is longitude and is between 0 and 2*pi.  The north pole has colatitude 0 and arbitrary longitude (although for most cases, use phi=), and longitudes are measured in anticlockwise direction.  Somew functions allow data to be given using Cartesian coordinates (that is, (x,y,z)).
}
\details{
\tabular{ll}{
Package: \tab spherstat\cr
Type: \tab Package\cr
Version: \tab 1.7-2\cr
Date: \tab 2016-02-02\cr
License: \tab GPL (>= 2.0) \cr
LazyLoad: \tab yes\cr
}
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
Adrian Baddeley
}

\note{
Many functions in these package are the analogues for point processes on the sphere of functions available in \pkg{spatstat} for point processes in \eqn{R^2}{R^2} or \eqn{R^3}{R^3}.  Hence elements of the code and help pages in this package have been taken from corresponding material in \pkg{spatstat}, with the permission of A. J. Baddeley.  This enables the code in \pkg{spatstat} to be both efficient and give similar output to corresponding functions in \pkg{spatstat}.  It is hoped that this will minimise or remove any confusion for users of both \pkg{spatstat} and \pkg{spherstat}.
}

\keyword{"spatial statistics"}

\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.Gest <- Gsphere(X=sph.pp, win=sph)
}
