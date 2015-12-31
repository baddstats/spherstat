\name{intensitysph}
\alias{intensitysph}
\title{
  Homogeneous Intensity of a Dataset
}
\description{
  A function for computing the intensity of a spatial dataset or spatial
  point process model, under the assumption of homogeneity. 
}
\usage{
 intensitysph(X, win, ...)
}
\arguments{
  \item{X}{
    The locations of the points in the point pattern.  An object of
    class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 column matrix
    with locations given in spherical coordinates, or a 3 column matrix
    with locations given in Cartesian coordinates. 
  }
  \item{win}{
    The window.  An object of type \code{\link{sphwin}}.  Only required
    if \code{X} is not an object of class \code{\link{sp2}} or
    \code{\link{sp3}}. 
  }
  \item{\dots}{
    Extra information.  Not used now, but probably needed when I expand
    this to other point processes 
  }
}
\details{
  Estimates the intensity by dividing the number of points in \code{X}
  by the area of \code{win}. 
}
\value{
A numeric.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{intensity}} in \pkg{spatstat}, which is
  the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements
  of this help page have been taken from \code{\link[spatstat]{intensity}}
  with the 
  permission of A. J. Baddeley.  This enables the information on this help
  page to be consistent with that for \code{\link[spatstat]{intensity}}.  It is
  hoped that this will minimise or remove any confusion for users of both
  \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link{rpoispp.sphwin}}
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sphwin(win=sph, lambda=10)
intensitysph(X=rps, win=sph)
}
\keyword{spatial}
\keyword{nonparametric}
