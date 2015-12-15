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
convert.globe(X)
convert.sphereplot(X)
}
\arguments{
  \item{X}{
    A 2 column matrix containing spherical coordinates of points to be plotted
  }
}
\value{
  A 2 column matrix that can be plotted using functions in the package
  \pkg{globe} (if \code{convert.globe} used) or \pkg{sphereplot} (if
  \code{convert.sphereplot} used) 
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  \code{convert.sphereplot} does not need the radius of the sphere; you
  specify the radius when making calls to \code{rgl.sphgrid} and
  \code{rgl.sphpoints} that plot the data. 
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
}
# if(require(globe)) {
# X <- rpoispp.sphwin(win=sphere, lambda=150, rad=1)
# X.globe <- convert.globe(X)
# globeearth(eye=northpole, top=nedlands)
# globepoints(loc=X.globe, eye=northpole, top=nedlands)
# }
}
\keyword{spatial}
\keyword{manip}
