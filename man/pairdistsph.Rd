\name{pairdist.sph}
\alias{pairdistsph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise distances
}
\description{
Computes the matrix of distances between all pairs of 'things' in a dataset, where the 'things' are on (a subset of) a sphere.
}
\usage{
pairdistsph(X, rad = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Locations of points.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 column matrix, with locations in spherical coordinates, or a 3 column matrix, with locations in Cartesian coordinates.
}
  \item{rad}{
Optional.  The radius of the sphere, default value 1.  Not required if \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}.
}
}
\details{
Given a dataset (representing a point pattern on a sphere, or a subset of a sphere) \code{pairdistsph} computes the minimum great circle distance between each pair of 'points' in the dataset, and returns a matrix containing these distances.

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius will be used rather than that specified in \code{rad}.
}
\value{
A square matrix whose \code{[i,j]} entry is the distance between the 'things' numbered i and j. 
}

\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{pairdist}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{pairdist}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{pairdist}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{nndistsph}} (nearest neighbour distances), \code{\link{gcdist}} (great circle distance between two points), \code{\link{bdist.sphwin}} (boundary distances)
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sp2(win=sph, lambda=1)
pairdistsph(points=rps)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
