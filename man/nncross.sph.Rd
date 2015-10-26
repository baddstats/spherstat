\name{nncross.sph}
\alias{nncrosssph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearest neighbour distances between 2 patterns.
}
\description{
Given two point patterns X and Y on a (subset of) a sphere, finds the nearest neighbour in Y of each point of X.
}
\usage{
nncrosssph(X, Y, rad=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X, Y}{
Locations of points.  Objects of class \code{\link{sp2}} or \code{\link{sp3}}; 2 column matrices, with locations in spherical coordinates; or 3 column matrices, with locations in Cartesian coordinates.  \code{X} and \code{Y} need not be given in the same format.
}

\item{rad}{
The radius of the sphere.
}
}
\details{
Given two point patterns X and Y this function finds, for each point of X, the nearest point of Y. This is done using great circle distances. 

The return value is a numeric, the ith value of which is the great circle distance from ith point in \code{X} to the nearest point in \code{Y}.

Note that this function is not symmetric in X and Y. To find the nearest neighbour in X of each point in Y, where Y is a point pattern, use nncross(Y,X). 

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation.  If \code{X} is a matrix, and \code{Y}  is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation.  If both \code{X} and \code{Y} are objects of class \code{\link{sp2}} or \code{\link{sp3}}, they must have the same radius or an error will be generated.
}
\value{
Numeric containing the nearest neighbour distances for each point. 
}

\author{
Tom Lawrence tjlawrence@bigpond.com
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{link{nndistsph}} (nearest neighbour distances for a single sample), \code{\link{crossdistsph}} (pairwise distances between two samples), \code{\link{pairdistsph}} (pairwise distances for a single sample), \code{\link{gcdist}} (great circle distance between two points), \code{\link{bdist.sphwin}}(boundary distances)
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sphere(win=sph, lambda=1)
rps2 <- rpoispp.sphere(win=sph, lambda=1)
X <- sp2(X=rps, win=sph)
Y <- sp2(X=rps2, win=sph)
nncrosssph(X=X, Y=Y)

nncrosssph(X=rps, Y=rps2, rad=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
