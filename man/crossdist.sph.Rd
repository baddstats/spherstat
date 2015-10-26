\name{crossdist.sph}
\alias{crossdistsph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise distances
}
\description{
Computes the distances between pairs of 'things' taken from two different datasets where the 'things' are on (a subset of) a sphere.
}
\usage{
crossdistsph(X, Y, rad = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X, Y}{
Locations of points.  Objects of class \code{\link{sp2}} or \code{\link{sp3}}; 2 column matrices, with locations in spherical coordinates; or 3 column matrices, with locations in Cartesian coordinates.  \code{X} and \code{Y} need not be given in the same format.
}
  \item{rad}{
The radius of the sphere, default value 1.  Only required when neither \code{X} nor \code{Y} are objects of class \code{\link{sp2}} or \code{\link{sp3}}.
}
}
\details{
Given two datasets \code{X} and \code{Y} (representing two point patterns on the same (subset of a) sphere) \code{crossdist.sph} computes the great circle distance from each thing in the first dataset to each thing in the second dataset, and returns a matrix containing these distances.

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation.  If \code{X} is a matrix, and \code{Y}  is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation.  If both \code{X} and \code{Y} are objects of class \code{\link{sp2}} or \code{\link{sp3}}, they must have the same radius or an error will be generated.
 
}
\value{
A matrix whose \code{[i,j]} entry is the distance from the i-th thing in the first dataset to the j-th thing in the second dataset. 
}

\author{
Tom Lawrence tjlawrence@bigpond.com
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{nncrosssph}} (nearest neighbour distances between two samples), \code{\link{nndistsph}} (nearest neighbour distances), \code{\link{gcdist}} (great circle distance between two points), \code{\link{bdist.sphwin}} (boundary distances)
}
\examples{
sph <- sphwin(type="sphere")
rps1 <- rpoispp.sp2(win=sph, lambda=1)
rps2 <- rpoispp.sp2(win=sph, lambda=1)
crossdistsph(X=rps1, Y=rps2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
