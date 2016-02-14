\name{convert.globe}
\alias{convert.globe}
\alias{convert.sphereplot}
\title{
Convert coordinates for plotting
}
\description{
Convert a matrix of spherical coordinates into a format that can be
plotted using functions available in either the \pkg{globe} or \pkg{sphereplot}
packages 
}
\usage{
convert.globe(X, inverse=FALSE, sp.dim=2)
convert.sphereplot(X, inverse=FALSE, sp.dim=2)
}
\arguments{
  \item{X}{
    If \code{inverse=TRUE}, a 2 or 3 column matrix containing spherical coordinates of points to be plotted.  
    If \code{inverse=FALSE} and using \code{convert.globe}, a list of length 2 with values \code{lat} and \code{lon}.
    If \code{inverse=FALSE} and using \code{convert.sphereplot}, a 2 column matrix or length 2 numeric.
  }
  \item{inverse}{
   Logical.  If TRUE, then conversion is to coordinates that can be used in \pkg{spherstat}.  If FALSE, then converstion is to coordinates that  can be used in \pkg{globe} (\code{convert.globe}) or \pkg{sphereplot} (code{convert.sphereplot}).
  }
  \item{sp.dim}{
   If \code{inverse=TRUE}, takes value 2 or 3 to respectively produce a 2 or 3 column matrix.  Ignored if \code{inverse=FALSE}.
  }
}
\value{
  If \code{inverse=FALSE}, a 2 column matrix that can be plotted using functions in the package
  \pkg{globe} (if \code{convert.globe} used) or \pkg{sphereplot} (if
,  \code{convert.sphereplot} used).  If \code{inverse=TRUE}, a 2 column matrix of locations in spherical coordinates (if \code{sp.dim=2}) or a 3 column matrix of locations in Cartesian coordinates (if \code{sp.dim=3}).
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  \code{convert.sphereplot} does not need the radius of the sphere; you
  specify the radius when making calls to \code{rgl.sphgrid} and
  \code{rgl.sphpoints} that plot the data.
   Also, if coordinates are converted using either function, and then converted back again with the same function, the doubly converted coordinates may differ from the orignal coordinates by a minute quantity (i.e. in the vicinity of 10^-16).
}
\seealso{
  \code{\link{convert2}}, \code{\link{convert3}}
}
\examples{
sphere <- sphwin(type="sphere")
X <- rpoispp.sphwin(win=sphere, lambda=150)
X.sphereplot <- convert.sphereplot(X)
if(require(sphereplot)) {
  rgl.sphgrid(longtype="D")
  rgl.sphpoints(long=X.sphereplot, radius=1)
Xback <- convert.sphereplot(X.sphereplot, inverse=TRUE, sp.dim=2)
}
if(require(globe)) {
X <- rpoispp.sphwin(win=sphere, lambda=150, rad=1)
X.globe <- convert.globe(X)
globeearth(eye=place("northpole"), top=place("nedlands"))
globepoints(loc=X.globe, eye=place("northpole"), top=place("nedlands"))
Xback <- convert.globe(X.globe, inverse=TRUE, sp.dim=3)
}
}
\keyword{spatial}
\keyword{manip}
