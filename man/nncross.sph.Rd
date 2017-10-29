\name{nncross.sph}
\alias{nncrosssph}
\title{
  Nearest Neighbour Distances Between Two Point Patterns on the Sphere
}
\description{
  Given two point patterns \code{X} and \code{Y} on a (subset of) a
  sphere, finds the nearest neighbour in \code{Y} of each point of \code{X}. 
}
\usage{
nncrosssph(X, Y, rad=1)
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
    The radius of the sphere.
  }
}
\details{
  Given two point patterns \code{X} and \code{Y} this function finds,
  for each point of \code{X}, the nearest point of \code{Y}. This is
  done using great circle distances.  

  The return value is a numeric, the \code{i}th value of which is the
  great circle distance from the \code{i}th point in \code{X} to the
  nearest point in \code{Y}. 

  Note that this function is not symmetric in \code{X} and \code{Y}. To
  find the nearest neighbour in \code{X} of each point in \code{Y},
  where \code{Y} is a point pattern, use \code{nncross.sph(Y,X)}.  

  If \code{X} is an object of class \code{\link{sp2}} or
  \code{\link{sp3}}, then its radius is used in this calculation.  If
  \code{X} is a matrix, and \code{Y}  is an object of class
  \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in
  this calculation.  If both \code{X} and \code{Y} are objects of class
  \code{\link{sp2}} or \code{\link{sp3}}, they must have the same radius
  or an error will be generated. 
}
\value{
  Numeric containing the nearest neighbour distances for each point. 
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\seealso{
  \code{\link{nndistsph}} (nearest neighbour distances for a single
  sample), \code{\link{crossdistsph}} (pairwise distances between two
  samples), \code{\link{pairdistsph}} (pairwise distances for a single
  sample), \code{\link{gcdist}} (great circle distance between two
  points), \code{\link{bdist.sphwin}}(boundary distances), \code{\link[spatstat]{nncross}}
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{nncross}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
  elements of the code for \code{nncross.sph} and help page have been taken
  from \code{\link[spatstat]{nncross}} with the permission of
  A. J. Baddeley.  This enables the code to be highly efficient and give
  corresponding output to, and for the information on this help page to
  be consistent with that for the function \code{\link[spatstat]{nncross}}.
  It is hoped that this will minimise or remove any confusion for users
  of both \pkg{spatstat} and \pkg{spherstat}. 
}

\examples{
sph <- sphwin(type="sphere")
X <- rpoispp.sphwin(win=sph, lambda=1)
Y <- rpoispp.sphwin(win=sph, lambda=1)
nncrosssph(X, Y)
Xmat <- rpoispp.sphwin(win=sph, lambda=1, as.sp=FALSE)
Ymat <- rpoispp.sphwin(win=sph, lambda=1, as.sp=FALSE)
nncrosssph(Xmat, Ymat, rad=1)
}
\keyword{spatial}
\keyword{math}

