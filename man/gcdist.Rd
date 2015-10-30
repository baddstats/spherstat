\name{gcdist}
\alias{gcdist}
\alias{gcdistPaired}
\title{
  Great Circle Distances
}
\description{
  Calculate the minimum great circle distance between points on a sphere
}
\usage{
gcdist(x, y, rad = 1)
gcdistPaired(x, y, rad = 1)
}
\arguments{
  \item{x,y}{
    The locations of the two points for which the distance is to be
    calculated.  Each must be a matrix with 2 columns (with location in
    spherical coordinates) or or 3 (location in Cartesian coordinates),
    or objects of class \code{\link{sp2}} or \code{\link{sp3}}.  Note
    \code{x} and \code{y} need not be in the same format. 
  }
  \item{rad}{
    Optional.  The radius of the sphere, default is 1.
  }
}
\details{
  For any two points on a sphere, the shortest path (also on the sphere)
  is the minor arc of a great circle (that is, a circle that is concentric
  with the sphere, and is a subset of the sphere).

  The function \code{gcdist} calculates the great circle distance
  between each pair of points in \code{x} and \code{y}, returning a
  matrix of size \code{m * n} where \code{m} is the number of
  points in \code{x} and \code{n} is the number of points in \code{y}.

  The function \code{gcdistPaired} calculates the great circle distance
  between corresponding points in \code{x} and \code{y} (that is,
  between the \code{i}th point in \code{x} and the \code{i}th point in
  \code{y} for each index \code{i}), returning a
  numeric vector. The point patterns \code{x} and \code{y} should
  contain equal numbers of points.
}
\value{
  A numeric matrix or vector
  giving the minimum great circle distance between the given pairs of points.
}
\author{
Tom Lawrence
}
\seealso{
\code{\link{nndistsph}} (nearest neighbour distances),
\code{\link{pairdistsph}} (pairwise distances),
\code{\link{bdist.sphwin}} (boundary distances) 
}
\examples{
x1 <- matrix(c(pi/2, 0), ncol=2)
y1 <- matrix(c(pi/2, pi/4), ncol=2)
gcdist(x=x1, y=y1)
}
\keyword{spatial}
\keyword{math}
