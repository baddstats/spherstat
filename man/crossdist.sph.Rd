\name{crossdist.sph}
\alias{crossdistsph}
\title{
  Pairwise Distances
}
\description{
Computes the distances between pairs of points on the sphere, taken from two
different datasets
}
\usage{
crossdistsph(X, Y, rad = 1)
}
\arguments{
  \item{X, Y}{
    Locations of points.  Objects of class \code{\link{sp2}} or
    \code{\link{sp3}}; 2 column matrices, with locations in spherical
    coordinates; or 3 column matrices, with locations in Cartesian
    coordinates.  \code{X} and \code{Y} need not be given in the same
    format. 
  }
  \item{rad}{
    The radius of the sphere, default value 1.  Only required when
    neither \code{X} nor \code{Y} are objects of class \code{\link{sp2}}
    or \code{\link{sp3}}. 
  }
}
\details{
  Given two datasets \code{X} and \code{Y} (representing two point
  patterns on the same (subset of a) sphere) \code{crossdist.sph} computes
  the great circle distance from each point in the first dataset to each
  point in the second dataset, and returns a matrix containing these
  distances. 

  If \code{X} is an object of class \code{\link{sp2}} or
  \code{\link{sp3}}, then its radius is used in this calculation.  If
  \code{X} is a matrix, and \code{Y}  is an object of class
  \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in
  this calculation.  If both \code{X} and \code{Y} are objects of class
  \code{\link{sp2}} or \code{\link{sp3}}, they must have the same radius
  or an error will be generated. 
}
\value{
  A matrix whose \code{[i,j]} entry is the distance from the \code{i}-th
  point in the first dataset to the \code{j}-th point in the second dataset. 
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat.geom]{crossdist}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
  elements of the code for \code{crossdist.sph} and help page have been taken
  from \code{\link[spatstat.geom]{crossdist}} with the permission of
  A. J. Baddeley.  This enables the code to be highly efficient and give
  corresponding output to, and for the information on this help page to
  be consistent with that for the function \code{\link[spatstat.geom]{crossdist}}.
  It is hoped that this will minimise or remove any confusion for users
  of both \pkg{spatstat} and \pkg{spherstat}. 
}

\seealso{
  \code{\link{nncrosssph}} (nearest neighbour distances between two
  samples), \code{\link{nndistsph}} (nearest neighbour distances),
  \code{\link{gcdist}} (great circle distance between two points),
  \code{\link{bdist.sphwin}} (boundary distances), \code{\link[spatstat.geom]{crossdist}}
}
\examples{
sph <- sphwin(type="sphere")
rps1 <- rpoispp.sphwin(win=sph, lambda=1)
rps2 <- rpoispp.sphwin(win=sph, lambda=1)
crossdistsph(X=rps1, Y=rps2)
}
\keyword{spatial}
\keyword{math}

