\name{sppatt}
\alias{sp2}
\alias{sp3}
\title{
  Create a Point Pattern on a (Subset of) a Sphere
}
\description{
  Creates an object of class \code{"sp2"} or \code{"sp3"} representing a
  point pattern on a (subset of) a sphere  
}
\usage{
sp2(X, win=sphwin(type="sphere"), check=TRUE)
sp3(X, win=sphwin(type="sphere"), check=TRUE)
}
\arguments{
  \item{X}{
    The locations of points in the pattern.  \code{sp2} requires a 2
    column matrix, with locations given in spherical coordinates (see
    \code{\link{spherstat-package}}), \code{sp3} requires a 3 column
    matrix, with locations given in Cartesian coordinates.  See Details
    for further information on how to specify locations. 
  }
  \item{win}{
    The observation window, an object of class \code{\link{sphwin}}.
  }
   \item{check}{
     If \code{TRUE} (which is the default setting),
     checks that all points in \code{X} are in the window defined by
     \code{win}.  If FALSE, no checks are performed.  See Details.  
   }
}
\details{
  An object of class \code{sp2} or \code{sp3} is the most efficient way
  of giving information about the locations of points on a sphere, and
  the window in which they are observed.  Functions in this package
  accept objects of these types as input (even when one of the elements
  e.g. the locations of points, or the window, is irrelevant), they can
  instead accept the relevant information provided in the same way as
  specified above e.g. if \code{win} is specified as an object of class
  \code{"sphwin"}.  This being said, \code{sp2} and \code{sp3} are
  strongly recommended as the classes of objects given as input in this
  package because their creation includes a check to ensure all
  specified points are in \code{win} and in the required format. 

  For both \code{sp2} and \code{sp3}, the locations of points (\code{X})
  must be given in matrix form, with each row specifying the location of
  a different point.  \code{sp2} requires a 2 column matrix with
  locations given in spherical coordinates; the first column being
  colatitudes, taking values in the range \eqn{[0, \pi]}{[0, pi]},
  where \eqn{0} is the `north pole'; the second column being longitudes,
  taking values in the range \eqn{[0, 2\pi)}{[0, 2*pi)}, where
  \eqn{0} is the `Greenwich meridian', and the
  longitude is measured in an anticlockwise direction from this point. 

  \code{sp3} requires a three column matrix; the first, second and third
  columns specifying the x, y and z coordinates of points.  A check is
  automatically performed to ensure that distance between the centre of
  the sphere (assumed to be (0,0,0)), the point (x,y,z) is within
  \code{10^(-11)} units of the surface of the sphere. 

  By default, when objects of class \code{sp2} and \code{sp3} are
  created, a check is performed to ensure that all points listed in
  \code{X} are within the observation window defined by \code{win}.  The
  argument \code{check} allows the user to determine whether or not such
  checks occur, and was introduced to reduce the run-time of some other
  functions in this package.  It is very strongly recommended that this
  argument always be set to \code{TRUE} (the default setting) since if
  an object of class \code{sp2} or \code{sp3} is created such that some
  points are outside the observation window, and that object is passed
  to other functions then either errors will be generated, or the output
  is likely to be invalid. 
}
\value{
  A one sentence description of the window.
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\seealso{
  \code{\link{spherstat-package}} for detail on the spherical coordinate
  system used in this package 
}
\examples{
X <- matrix(c(1,0,0, 0,0,1, 0,1,0), nrow=3, ncol=3, byrow=TRUE)
win <- sphwin(type="sphere")
sp3(X=X, win=win)

Y <- matrix(c(0,0, pi/2, 0, pi/2, pi, pi, 0), nrow=4, ncol=2, byrow=TRUE)
sp2(X=Y, win=win)

}
\keyword{spatial}
\keyword{manip}
