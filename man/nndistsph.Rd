\name{nndistsph}
\alias{nndistsph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearest neighbour distances
}
\description{
Computes the distance from each point to its nearest neighbour in a point pattern on (a subset of) a sphere.
}
\usage{
nndistsph(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Locations of points.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, a 2 column matrix, with locations in spherical coordinates, or a 3 column matrix, with locations in Cartesian coordinates.
}
}
\details{
This function computes the minimum great circle distance from each point in a point pattern on a (subset of) a sphere to its nearest neighbour (the nearest other point of the pattern).

To use the nearest neighbour distances for statistical inference, it is often advisable to use the edge-corrected empirical distribution, computed by \code{\link{Gsphere}}. 

To find the nearest neighbour distances from one point pattern to another point pattern, use nncross.

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius will be used rather than that specified in \code{rad}.
}
\value{
Numeric vector or matrix containing the nearest neighbour distances for each point. 
}

\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{nndist}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{nndist}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{nndist}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pairdistsph}} (pairwise distances for a single sample), \code{\link{gcdist}} (great circle distance between two points), \code{\link{bdist.sphwin}}(boundary distances)
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sp2(win=sph, lambda=1)
nndistsph(X=rps)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
