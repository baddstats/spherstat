\name{summary.sphwin}
\alias{summary.sphwin}
\title{
Summary of a Spatial Window that is a (Subset of) a Sphere
}
\description{
Prints a useful description of a window object. 
}
\usage{
\method{summary}{sphwin}(object, ...)
}
\arguments{
  \item{object}{ 
    Window (object of class "\code{sphwin}").
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
A useful description of the window \code{object} is printed. 

This is a method for the generic function \code{\link{summary}}. 
}

\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for windows on the sphere of the
  function \code{\link[spatstat.geom]{summary.owin}} in \pkg{spatstat}, which
  is the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
  elements of this help page have been taken from
  \code{\link[spatstat.geom]{summary.owin}} with the permission of
  A. J. Baddeley.  This enables the information on this help page to be
  consistent with that for \code{\link[spatstat.geom]{summary.owin}}.  It is
  hoped that this will minimise or remove any confusion for users of
  both \pkg{spatstat} and \pkg{spherstat}. 
}

\seealso{
\code{\link{summary.sp2}}, \code{\link{summary.sp3}},
\code{\link{print.sp2}}, \code{\link{print.sp3}},
\code{\link{print.sphwin}} 
}
\examples{
sph <- sphwin(type="sphere")
print(sph)
summary(sph)

sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
print(sphcap)
summary(sphcap)
}
\keyword{spatial}
\keyword{methods}

