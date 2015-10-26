\name{convertspherCart}
\alias{convert2}
\alias{convert3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between spherical and Cartesian coordinates
}
\description{
convert2 converts from Cartesian to spherical coordinates, whereas convert3 converts from spherical to Cartesian coordinates.
}
\usage{
convert2(points, rad = 1)
convert3(points, rad = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
The coordinates of the points to be converted.  See Details.  
}
  \item{rad}{
The radius of the sphere, by default set to 1.  Not required if \code{points} is an object of class \code{\link{sp2}} or \code{\link{sp3}}.
}
}
\details{
For \code{convert2}, \code{points} must contain locations of a point or points in Cartesian coordinates.  This can be an object of class \code{\link{sp3}} or a length 3 vector or a 3 column matrix.

For \code{convert3}, \code{points} must contain locations of a point or points in spherical coordinates.  This can be an object of class \code{\link{sp2}} or a length 2 vector or a 2 column matrix.
}
\value{
The output of convert2 is either an object of class \code{\link{sp2}} (if \code{points} is an object of class \code{\link{sp2}}), or a 2 column matrix giving the locations of points in spherical coordinates.  The output of convert2 is either an object of class \code{\link{sp2}} (if \code{points} is an object of class \code{\link{sp3}}), or a 3 column matrix giving the locations of points in Cartesian coordinates.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{convert.globe}}, \code{\link{convert.sphereplot}}
}
\examples{
sphere <- sphwin(type="sphere")
X <- rpoispp.sp2(sphere, lambda=20)
Y <- convert3(X)

i3 <- matrix(c(1,0,0,0,1,0,0,0,1), nrow=3, ncol=3, byrow=TRUE)
convert2(i3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
