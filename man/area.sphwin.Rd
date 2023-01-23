\name{area.sphwin}
\alias{area.sphwin}
\title{
Calculate the area on the surface of a sphere
}
\description{
Calculates the area of a sphere or a subset of the sphere (cap, band, wedge or polygon)
}
\usage{
 \method{area}{sphwin}(w)
}
\arguments{
  \item{w}{
    The window for which the area is calculated.
    An object of type \code{"sphwin"}.
  }
}
\value{
  The area of the sphere is returned as an object of type \code{numeric}.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat.geom]{area.owin}} in \pkg{spatstat}, which is
  the corresponding function for point processes in \eqn{R^2}. 
}
\seealso{
  \code{\link{eroded.areas.sphwin}}, \code{\link[spatstat.geom]{area.owin}}
}
\examples{
sph <- sphwin(type="sphere")
area.sphwin(w=sph)

sphcap <- sphwin(type="sphere", param=c(0, pi/3))
area.sphwin(w=sphcap)
}
\keyword{spatial}
\keyword{math}
\keyword{methods}
