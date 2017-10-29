\name{bdist.sphwin}
\alias{bdist.sphwin}
\title{
  Boundary Distances
}
\description{
  Calculate the minimum great circle distance from a point (or each of a set of
  points) to the boundary of a window on the sphere. 
}
\usage{
bdist.sphwin(X, win = sphwin(type = "sphere"))
}
\arguments{
  \item{X}{
    The locations of the points for which the boundary distance is
    calculated.  Either a 2 column matrix with locations of points in
    spherical coordinates, a 3 column matrix with locations of points in
    Cartesian coordinates, or an object of class \code{"sp2"} or
    \code{"sp3"}. 
  }
  \item{win}{
    The window.  An object of class \code{"sphwin"}.
  }
}
\details{
  This function uses the dot product to calculates\ the minimum great circle distance from a point (or each of
  a set of points) to the boundary of a region on the sphere that can be
  defined using \code{\link{sphwin}}. 
}
\value{
  A numeric of length equal to the number of points in \code{X}.  The
  first value is the minimum boundary distance for the point in the first
  row of \code{X}, and so on. 
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point patterns on the sphere of the
function \code{\link[spatstat]{bdist.points}} in \pkg{spatstat}. 
}
\seealso{
  \code{\link[spatstat]{bdist.points}}, \code{\link{nndistsph}} (nearest
  neighbour distances), \code{\link{pairdistsph}} (pairwise distances),
  \code{\link{gcdist}} (great circle distance between two points) 
}
\examples{
sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
rpp.sc <- rpoispp.sphwin(win=sphcap, lambda=10)
bdist.sphwin(X=rpp.sc, win=sphcap)
}
\keyword{spatial}
\keyword{math}
