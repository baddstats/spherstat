\name{pairdist.sph}
\alias{pairdistsph}
\title{
  Pairwise Distances Between Points on a Sphere
}
\description{
  Computes the matrix of distances between all pairs of points on a sphere.
}
\usage{
pairdistsph(X, rad = 1)
}
\arguments{
  \item{X}{
    Locations of points.  An object of class \code{\link{sp2}} or
    \code{\link{sp3}}, or a 2 column matrix, with locations in spherical
    coordinates, or a 3 column matrix, with locations in Cartesian
    coordinates. 
  }
  \item{rad}{
    Optional.  The radius of the sphere, default value 1.  Not required
    if \code{X} is an object of class \code{\link{sp2}} or
    \code{\link{sp3}}. 
  }
}
\details{
  Given a dataset (representing a point pattern on a sphere, or a subset
  of a sphere) \code{pairdistsph} computes the minimum great circle
  distance between each pair of points in the dataset, and returns a
  matrix containing these distances. 

  If \code{X} is an object of class \code{\link{sp2}} or
  \code{\link{sp3}}, then its radius will be used rather than that
  specified in \code{rad}. 
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance between the
  points numbered i and j.  
}

\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\note{
This function is the analogue for point processes on the sphere of the
function \code{\link[spatstat.geom]{pairdist}} in \pkg{spatstat}, which is
the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements
of this help page have been taken from \code{\link[spatstat.geom]{pairdist}} with the
permission of A. J. Baddeley.  This enables the information on this help
page to be consistent with that for \code{\link[spatstat.geom]{pairdist}}.
It is hoped that this will minimise or remove any confusion for users of both
\pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
 \code{\link{nndistsph}} (nearest neighbour distances),
 \code{\link{gcdist}} (great circle distance between two points),
 \code{\link{bdist.sphwin}} (boundary distances), \code{\link[spatstat.geom]{pairdist}} 
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sphwin(win=sph, lambda=1)
pairdistsph(rps)
}
\keyword{spatial}
\keyword{math}

