\name{rotations}
\alias{rot.matrix}
\alias{rot.sphere}
\title{
  Rotating the Sphere
}
\description{
  Generate the matrix for a specified rotation of the sphere
  (\code{rot.matrix}) and determine the locations of a set of points under
  that rotation (\code{rot.sphere}). 
}
\usage{
rot.matrix(northpole, rad=1)
rot.sphere(points, northpole, rad = 1, inverse = FALSE)
}
\arguments{
  \item{points}{
    The coordinates of points to be rotated.  An object of class
    \code{\link{sp2}} or \code{\link{sp3}}, or  a length 2 numeric or 2
    column matrix (points given in spherical coordinates) or a length 3
    numeric or 3 column matrix (points given in Cartesian coordinates). 
  }
  \item{northpole}{
    The coordinates of a reference point, given in spherical coordinates
    as a length 2 numeric or 2 column, 1 row matrix. See argument
    \code{inverse} for more information. 
  }
  \item{rad}{
    Optional.  The radius of the sphere, by default set to 1.  Not required
    in \code{rot.sphere} if \code{points} is an object of class
    \code{\link{sp2}} or \code{\link{sp3}}. 
  }
  \item{inverse}{
    Logical.  Takes value \code{FALSE} for the rotation that sends the point
    given in \code{northpole} to the north pole (i.e. \code{c(0,0)}); or
    \code{TRUE} for the inverse of that rotation. 
  }
}
\details{
  It is sometimes necessary to find the locations of points on a sphere
  under a specified rotation; for example when one of the given points
  is rotated to the north pole.  The benefit of this is that some
  spherical geometric and trigonometric functions can be more easily
  evaluated, which is useful in the context of spatial point process
  analysis. 

  The function \code{rot.sphere} invokes \code{rot.matrix} to define the
  required rotation matrix; this rotation sends the point specified in
  \code{northpole} to the north pole (if \code{inverse=FALSE}) along the
  shortest great circle arc, or vice versa (if \code{inverse=TRUE}).
  However, this rotation may not be entirely suitable for all uses since
  the longitude at colatitudes \code{0, pi} are arbitrary, and are
  usually set to 0 for the purposes of calculations, including in the
  standard matrix rotations that \code{rot.matrix} generates matrices
  for, and \code{rot.sphere} performs.  As such, the rotation that is
  performed by \code{rot.sphere} may not be the desired rotation:
  e.g. for similar reasons to that described in the previous paragraph,
  we may need to rotate the sphere such that a point \eqn{x_1}{x1} be rotated to
  the north pole, and a second point \eqn{x_2}{x2} to a semicircle of longitude
  \eqn{\phi_0}{phi0}, however the rotation that sends \eqn{x_1}{x1} to
  the north pole sends \eqn{x_2}{x2} to a different longitude,
  \eqn{\phi_1}{phi1}.  This can be easily overcome by 
  first performing the rotation that sends \eqn{x_1}{x1} to the north pole, then
  performing a second roation on the result, where
  \code{northpole=c(0, phi1-phi0)} for \code{phi1 > phi0},
  or \code{northpole=c(0, 2*pi - (phi1-phi0)} for \code{phi1 < phi0}. 

  Under this second rotation (i.e. where the colatitude of the point given
  in \code{northpole} has colatitude 0), the axis of rotation passes
  through the north and south poles as well as the origin, thus \eqn{x_1}{x1} is
  invariant under the rotation.  Note that a rotation where the colatitude
  of the point given in \code{northpole} is pi does not do the same, it
  sends the south pole to the north pole along the semicircle of longitude
  specified in \code{northpole}.  Thus for rotations under which the north
  and south poles are invariant, ensure that the colatitude of the point
  passed to \code{northpole} is 0. 

  For \code{rot.sphere}, if \code{points} is an object of class
  \code{\link{sp2}} or \code{\link{sp3}}, then its radius will be used
  rather than that specified in \code{rad}. 
}
\value{
  The output of \code{rot.matrix} is the matrix for the specified
  rotation of the sphere; it has 3 columns and 3 rows.
  The output of \code{rot.sphere} is a 2 column matrix, giving the
  locations of the points under the specified rotation, in spherical
  coordinates. 
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\examples{
rot.matrix(northpole=c(pi/3, pi/3), rad=1)

rot.sphere(points=c(pi/3, pi/3), northpole=c(pi/3,pi/3), rad=1, inverse=FALSE)
rot.sphere(points=c(pi/3, pi/3), northpole=c(pi/3,pi/3), rad=1, inverse=TRUE)
}
\keyword{spatial}
\keyword{math}
