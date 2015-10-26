\name{bdist.sphwin}
\alias{bdist.sphwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boundary distances
}
\description{
Calculate the minimum great circle distance from each of a set of points to the boundary of a window on the sphere.
}
\usage{
bdist.sphwin(X, win = sphwin(type = "sphere"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The locations of the points for which the boundary distance is calculated.  Either a 2 column matrix with locations of points in spherical coordinates, a 3 column matrix with locations of points in Cartesian coordinates, or an object of class \code{\link{sp2}} or \code{\link{sp3}}.
}
  \item{win}{
The window.  An object of class{\link{sphwin}}
}

}
\details{
Calculates the minimum great circle distance from a point (or each of a set of points) to the boundary of a region on the sphere that can be defined using \code{\link{sphwin}}.
}
\value{
A numeric of length equal to the number of points in \code{X}.  The first value is the minimum boundary distance for the point in the first row of \code{X}, and so on.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}

\note{
This function is the analogue for point processes on the sphere of the function \code{\link{perimeter}} in \code{\link{spatstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{perimeter}}, \code{\link{nndistsph}} (nearest neighbour distances), \code{\link{pairdistsph}} (pairwise distances), \code{\link{gcdist}} (great circle distance between two points)
}
\examples{
sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
rpp.sc <- rpoispp.sp2(win=sphcap, lambda=10)
bdist.sphwin(X=rpp.sc, win=sphcap)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
