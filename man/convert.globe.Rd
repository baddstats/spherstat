\name{convertcoords}
\alias{convert.globe}
\alias{convert.sphereplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert coordinates for plotting
}
\description{
Convert a matrix spherical coordinates into a format that can be plotted using functions available in either the globe or sphereplot packages
}
\usage{
convert.globe(X)
convert.sphereplot(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A 2 column matrix containing spherical coordinates of points to be plotted
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A 2 column matrix that can be plotted using functions in the package globe (if convert.globe used) or sphereplot (if convert.sphereplot used)
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
convert.sphereplot does not need the radius of the sphere; you specify the radius when making calls to rgl.sphgrid and rgl.sphpoints that plot the data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
convert2, convert3
}
\examples{
sphere <- sphwin(type="sphere")
X <- rpoispp.sp2(win=sphere, lambda=150)
X.sphereplot <- convert.sphereplot(X)
library(sphereplot)
rgl.sphgrid(longtype="D")
rgl.sphpoints(long=X.sphereplot, radius=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
% sphere <- sphwin(type="sphere")
% X <- rpoispp.sphere(win=sphere, lambda=150, rad=1)
% X.globe <- convert.globe(X)
% globeearth(eye=northpole, top=nedlands)
% globepoints(loc=X.globe, eye=northpole, top=nedlands)