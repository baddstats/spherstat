\name{in.W}
\alias{in.W}
\title{
  Determining If a Point is in a Region
}
\description{
  This function determines whether a given point or points is/are inside
  the given window, which is a (subset) of the sphere. 
}
\usage{
in.W(points, win)
}
\arguments{
  \item{points}{
    A 2 column matrix giving the locations of all points to be tested.
  }
  \item{win}{
    The window.  An object of type \code{\link{sphwin}}.
  }
}
\details{
  This function establishes whether each point in \code{points} is in
  the window specified in \code{win}.  A point is considered to be in
  the window if it is within or on the boundaries of the window,
  including the vertices.

  The algorithm used in the case that \code{win} is a spherical polygon is that developed by Bevis and Chatelain (1989).
}

\value{
  A logical of same length as the number of points being tested.  The
  first element is TRUE if the point whose coordinates are in the first
  row of \code{points} is in the window specified in \code{win}, and
  FALSE otherwise. 
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}

\section{Warning}{
  On any sphere, the boundaries that define a window \eqn{W} on the
  sphere also define a second shape \eqn{W_2}{W2}, that containing all the
  boundaries and vertices of \eqn{W}, and the complement of \eqn{W}.  Since
  the two shapes have boundaries and vertices in common, the set of all
  points that this function would identify as being outside \eqn{W} is
  not the set of all points in \eqn{W_2}{W2}.
} 

\references{
Bevis, M. and Chatelain J.-L. (1989) Locating a point on a spherical surface relative to a spherical polygon of arbitrary shape. \emph{Mathematical Geology} \bold{21}, no. 8, pp. 811--827.
}

\examples{
sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
sphband <- sphwin(type="band", param=c(pi/6, pi/2), ref=c(0,0))
rpp.band <- rpoispp.sphwin(win=sphband, lambda=10)
in.W(points=rpp.band, win=sphcap)
}
\keyword{spatial}
\keyword{math}

