\name{eroded.areas.sphwin}
\alias{eroded.areas.sphwin}
\title{
  Area of morphological erosions of a sphere or a subset of a sphere
}
\description{
  Computes the areas of successive morphological erosions of a window,
  where the window is a sphere or subset of a sphere. 
}
\usage{
  eroded.areas.sphwin(win = sphwin(type = "sphere"),
                      r=NULL, method = "exact", ...)
}
\arguments{
  \item{win}{
    The window.
  }
  \item{r}{
    Numeric vector of radii at which erosions will be performed.  Must have
    length 512, and must be specified if \code{win} is of type
    \code{quadrangle} or \code{polygon}.  If \code{r=NULL} and \code{win} is
    of type \code{band} or \code{wedge}, then the \code{r} will be generated
    such that \code{r} is between 0 and the smallest distance at which the
    area of the erosion is 0.  If \code{win} is of type \code{sphere}, then
    \code{r} is arbitrary (see Details) but for consistency with other
    functions is a length 512 numeric with \code{r} between 0 and \code{pi}
    times the radius of the sphere. 
  }
  \item{method}{
    (For spherical wedges only) specifies the method used to
    estimate/calculate the area morphological erosions of the window.
    Must be one of \code{"exact"}, \code{"integral"} or \code{"grid"} (see Details). 
  }
  \item{\dots}{
    (For spherical wedges where method=\code{"grid"}, and polygons only) Input to
    be passed to \code{\link{gridmat}}.  See Details below, then
    \code{\link{gridmat}}, for more information. 
  }
}
\details{
  This function computes the areas of the erosions of the window w by
  each of the radii \code{r[i]}.  

  The morphological erosion of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset consisting of all points \eqn{x \in W}{x in W}
  such that the distance
  from \eqn{x} to the boundary of \eqn{W} is greater than or equal to
  \eqn{r}. In other words it is the result of trimming a margin of width
  \eqn{r} off the set \eqn{W}.   

  The argument \eqn{r} should be a vector of nonnegative numbers. The argument
  w should be a window (an object of class \code{"sphwin"}, see
  \code{\link{sphwin}}).

  For the sphere, \code{eroded.areas.sph} returns a length 512 numeric
  where each value is the area of the sphere, since the erosions cannot
  exist because there are no boundaries on the sphere. 

  For the spherical wedge, it is recommended that the exact method of
  calculating the areas of the morphological erosions (method="exact")
  as it exactly calculates the areas and has the least computation time
  of the three methods. 

  Alternatively, there are two estimates of the morphological erosions
  of the spherical wedge.  The first estimator is "integral", which
  involves estimating the integral formula for the morphological
  erosion:
  \deqn{
    \int_{0}^{r} \mbox{length}(\partial (W_{-s})) ds
  }{
    integral[0, r] length(delta(W[-s])) ds
  }
  This is a relatively accurate estimator and has a similar computation
  time to the exact method.   

  The second estimator, which is also the estimator used for the
  spherical polygon, is \code{"grid"}, which involves superimposing a grid of
  cells of equal area and angular width; for each \eqn{r}, the total number
  of cells for which the centre is \eqn{\ge r}{>=r} from the boundary of the
  wedge is found and the result multipled by the area of a single cell.
  This method is the slowest and generates a curve which is only smooth
  when the grid is sufficiently fine.  This estimator requires further
  arguments to be sent to \code{\link{gridmat}}; for the spherical
  polygon, all arguments that \code{\link{gridmat}} needs are required
  i.e. \code{colats, lons, ncolat, nlon}.  For the spherical wedge, only
  \code{ncolat, nlon} are required; \code{eroded.areas.sph} calculates
  the values of \code{colats, lons}. 
}
\value{
  A vector containing the areas of morphological erosions of the sphere.
}
\author{
  Tom Lawrence  <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{eroded.areas}} in \pkg{spatstat}, which
  is the corresponding function for point processes in R^2.  Hence
  elements of this help page have been taken from
  \code{\link[spatstat]{eroded.areas}} with the permission of A. J
  Baddeley.  This enables the information on this help page to be
  consistent with that for \code{\link[spatstat]{eroded.areas}}.  It is
  hoped that this will minimise or remove any confusion for users of
  both \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link{area.sphwin}} (area of a sphere or subset of a sphere),
  \code{\link{bdist.sphwin}} (minimum distance from each of one or more
  points to the boundary of a window). 
}
\examples{
sph <- sphwin(type="sphere")
eroded.areas.sphwin(win=sph, r=seq(0, pi, length=512))

sphwedge <- sphwin(type="wedge", param=c(pi/3, 0), ref=c(0,0))
eroded.areas.sphwin(win=sphwedge, r=seq(0, pi/6, length=512), method="exact")
}
\keyword{spatial}
\keyword{math}
