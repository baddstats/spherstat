\name{sphsin}
\alias{sphsin}
\title{
  Spherical Sine Function
}
\description{
  Given the length of one side of a spherical triangle and size of the
  angle opposite, this function can calculate the length of another side
  of the triangle (given also the size of the angle opposite it), or the
  size of another angle in the triangle (given also the length of the
  side opposite it). 
}
\usage{
sphsin(d1, d2, theta1, theta2, rad = 1)
}
\arguments{
  \item{d1}{
    The length of the known side (opposite the angle of known size)
  }
  \item{d2}{
    If the size of an angle is to be calculated, the length of the side
    opposite that angle.  Otherwise, must be set to \code{NULL}. 
  }
  \item{theta1}{
    The size of the known angle (opposite the side of known length)
  }
  \item{theta2}{
    If the length of a side is to be calculated, the size of the angle
    opposite that side.  Otherwise, must be set to \code{NULL}. 
  }
  \item{rad}{
    The radius of the sphere
  }
}
\details{
  This function can be used to either calculate an angle in a spherical
  triangle when the length of the opposite side is knwon, the length of
  another side and size of the angle opposite it are all known.  This
  function can also be used to calculate the length of a side of a
  spherical triangle when the size of the angle opposite it, the size of
  another angle in the the triangle and length of the side opposite it
  are all known.  The calculation performed depends on whether \code{d2}
  or \code{theta2} is set to \code{NULL}; the variable set to
  \code{NULL} is the one for which this function calculates the value.
  If both \code{d2} and \code{theta2} (or neither) are set to
  \code{NULL}, an error is generated. 

  The distance(s) given to this function should be across the surface of
  the sphere (i.e. not Euclidean), the size of the angle (s) should be
  given in radians. 
}
\value{
  If \code{sphsin} is used to calculate the length of a side, that
  length is given across the surface of the sphere; if \code{sphsin} is
  used to calculate the size of an angle, it is given in radians. 
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\seealso{
  \code{\link{sphcos}} (spherical cosine function),
  \code{\link{sph.angles}} (calculate the sizes of all angles in a
  spherical polygon) 
}
\examples{
sphsin(d1=pi/2, d2=NULL, theta1=pi/2, theta2=pi/3, rad=1)

sphsin(d1=pi/2, d2=pi/3, theta1=pi/2, theta2=NULL, rad=1)}
}
\keyword{spatial}
\keyword{math}

