\name{sphcos}
\alias{sphcos}
\title{
The spherical cosine function
}
\description{
Use the spherical cosine function to calculate the length of a side or size of an angle in a spherical triangle.
}
\usage{
sphcos(d1, d2, d3, theta, rad = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d3}{
If the size of an angle is to be calculated, the length of the side opposite that angle.  Otherwise, must be set to \code{NULL}.
}
  \item{theta}{
If the length of a side is to be calculated, the size of the angle opposite that side.  Otherwise, must be set to \code{NULL}.
}
  \item{d1, d2}{
The lengths of the other two sides
}
  \item{rad}{
The radius of the sphere
}
}
\details{
This function can be used to either calculate an angle in a spherical triangle when the lengths of all three sides are known, or the length of a side of a spherical triangle when the size of the angle opposite it and the lengths of the other two sides are known.  The calculation performed depends on whether \code{d3} or \code{theta} is set to \code{NULL}; the variable set to \code{NULL} is the one for which this function calculates the value.  If both \code{d3} and \code{theta} (or neither) are set to \code{NULL}, an error is generated.

The distances given to this function should be across the surface of the sphere (i.e. not Euclidean), the size of the angle (if this function is used to calculate the length of the opposite side) should be given in radians.
}
\value{
If \code{sphcos} is used to calculate the length of a side, that length is given across the surface of the sphere; if \code{sphcos} is used to calculate the size of an angle, it is given in radians.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}

\seealso{
\code{\link{sphsin}} (spherical sine function), \code{\link{sph.angles}} (calculate the sizes of all angles in a spherical polygon)
}
\examples{
sphcos(d1=pi/2, d2=pi/2, d3=pi/3, theta=NULL, rad=1)

sphcos(d1=pi/2, d2=pi/2, d3=NULL, theta=pi/3, rad=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
