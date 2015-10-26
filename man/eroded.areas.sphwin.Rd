\name{eroded.areas.sphwin}
\alias{eroded.areas.sphwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Area of morphological erosions of a sphere or a subset of a sphere
}
\description{
Computes the areas of successive morphological erosions of a window, where the window is a sphere or subset of a sphere.
}
\usage{
eroded.areas.sphwin(win = sphwin(type = "sphere"), r=NULL, method = "exact", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{win}{
The window.
}
  \item{r}{
Numeric vector of radii at which erosions will be performed.  Must have length 512, and must be specified if \code{win} is of type \code{qaudrangle} or \code{polygon}.  If \code{r=NULL} and \code{win} is of type \code{band} or \code{wedge}, then the \code{r} will be generated such that \code{r} is between 0 and the smallest distance at which the area of the erosion is 0.  If \code{win} is of type \code{sphere}, then \code{r} is arbitrary (see Details) but for consistency with other functions is a length 512 numeric with \code{r} between 0 and \code{pi} times the radius of the sphere.
}
  \item{method}{
(For spherical wedges only) specifies the method used to estimate/calculate the area morphological erosions of the window.  Must be one of "exact", "integral" or "grid" (see Details).
}
  \item{\dots}{
(For spherical wedges where method="grid", and polygons only) Input to be passed to \code{\link{gridmat}}.  See Details below, then \code{\link{gridmat}}, for more information.
}
}
\details{
This function computes the areas of the erosions of the window w by each of the radii r[i]. 

The morphological erosion of a set W by a distance r > 0 is the subset consisting of all points x in W such that the distance from x to the boundary of W is greater than or equal to r. In other words it is the result of trimming a margin of width r off the set W. 

The argument r should be a vector of nonnegative numbers. The argument w should be a window (an object of class "sphwin", see help(sphwin.object). 

For the sphere, \code{eroded.areas.sph} returns a length 512 numeric where each value is the area of the sphere, since the erosions cannot exist because there are no boundaries on the sphere.

For the spherical wedge, it is recommended that the exact method of calculating the areas of the morphological erosions (method="exact") as it exactly calculates the areas and has the least computation time of the three methods.

Alternatively, there are two estimates of the morphological erosions of the spheriucl wedge.  The first estimator is "integral", which involves estimating the integral formula for the morphological erosion: 
\code{int_{0}^{r}length(delta (W_{-s}))ds}
This is a relatively accurate estimator and has a similar computation time to the exact method.  

The second estimator, which is also the estimator used for the spherical polygon, is "grid", which involves superimposing a grid of cells of equal area and angular width; for each $r$, the total number of cells for which the centre is \code{>=r} from the boundary of the wedge is found and the result multipled by the area of a single cell.  This method is the slowest and generates a curve which is only smooth when the grid is sufficiently fine.  This estimator requires further arguments to be sent to \code{\link{gridmat}}; for the spherical polygon, all arguments that \code{\link{gridmat}} needs are required i.e. \code{colats, lons, ncolat, nlon}.  For the spherical wedge, only \code{ncolat, nlon} are required; \code{eroded.areas.sph} calculates the values of \code{colats, lons}.


}
\value{
A vector containing the areas of morphological erosions of the sphere.
}
\author{
Tom Lawrence
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{eroded.areas}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{eroded.areas}} with the permission of A. J Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{eroded.areas}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
area.sphwin (area of a sphere or subset of a sphere), bdist.sphwin (minimum distance from each of one or more points to the boundary of a window).
}
\examples{
sph <- sphwin(type="sphere")
eroded.areas.sphqin(win=sph, r=seq(0, pi, length=512))

sphwedge <- sphwin(type="wedge", param=c(pi/3, 0), point=c(0,0))
eroded.areas.sphwin <- sphwin(win=sphwedge, r=seq(0, pi/6, length=512), method="exact")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
