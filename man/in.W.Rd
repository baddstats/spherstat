\name{in.W}
\alias{in.W}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determining if a point is in a region
}
\description{
This function determines whether a given point or points is/are inside the given window, which is a (subset) of the sphere.
}
\usage{
in.W(points, win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
A 2 column matrix giving the locations of all points to be tested.
}
  \item{win}{
The window.  An object of type \code{\link{sphwin}}.
}
}
\details{
This function establishes whether each point in \code{points} is in the window specified in \code{win}.  A point is considered to be in the window if it is within or on the boundaries of the window, including the vertices.
}

\value{
A logical of same length as the number of points being tested.  The first element is TRUE if the point whose coordinates are in the first row of \code{points} is in the window specified in \code{win}, and FALSE otherwise.
}

\author{
Tom Lawrence
}

\section{Warning }{On any sphere, the boundaries that define a window \code{W} on the sphere also define a second shape \code{W_2}, that containing all the boundaries and vertices of \code{W}, and the complement of W.  Since the two shapes have boundaries and vertices in common, the set of all points that this function would identify as being outside \code{W} is not the set of all points in \code{W_2}.}

\examples{
sphcap <- sphwin(type="band", param=c(0, pi/3), point=c(0,0))
sphband <- sphwin(type="band", param=c(pi/6, pi/2), point=c(0,0))
rpp.band <- rpoispp.sp2(win=sphband, lambda=10)
in.W(point=rpp.band, win=sphcap, rad=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
