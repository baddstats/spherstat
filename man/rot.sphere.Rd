\name{rotations}
\alias{rot.matrix}
\alias{rot.sphere}
\alias{rot.sphere2}
\title{
  Rotating the Sphere
}
\description{
  Generate the matrix for a specified rotation of the sphere
  (\code{rot.matrix}) and determine the locations of a set of points under
  that rotation (\code{rot.sphere}).  Also, rotate a set of points through a specified angle about a specified axis
}
\usage{
rot.matrix(northpole, rad=1)
rot.sphere(points, northpole, rad = 1, inverse = FALSE)
rot.sphere2(points, v, theta)
}
\arguments{
  \item{points}{
    The coordinates of points to be rotated.  An object of class
    \code{\link{sp2}} or \code{\link{sp3}}, or  a length 2 numeric or 2
    column matrix (points given in spherical coordinates) or a length 3
    numeric or 3 column matrix (points given in Cartesian coordinates). 
  }
  \item{northpole}{
    For \code{rot.matrix}, \code{rot.sphere}: the coordinates of a reference point, given in spherical coordinates
    as a length 2 numeric or 2 column, 1 row matrix. See argument
    \code{inverse} for more information. 
  }
  \item{rad}{
    Optional.  For \code{rot.matrix}, \code{rot.sphere}: the radius of the sphere, by default set to 1.  Not required
    in \code{rot.sphere} if \code{points} is an object of class
    \code{\link{sp2}} or \code{\link{sp3}}. 
  }
  \item{inverse}{
    Logical.  For \code{rot.matrix}, \code{rot.sphere}: takes value \code{FALSE} for the rotation that sends the point
    given in \code{northpole} to the north pole (i.e. \code{c(0,0)}); or
    \code{TRUE} for the inverse of that rotation. 
  }
  \item{v}{
    For \code{rot.sphere2}: the axis about which to rotate the points.  A matrix with 2 or 3 columns
  }
  \item{theta}{
   For \code{rot.sphere2}: the angle through which the points should be rotated, in a clockwise direction. A numeric of length 1, should be at least 0 but less than \code{2*pi}.
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
  northpole\eqn{-c(0, \phi_{1}-\phi_{0})}{=c(0, phi1-phi0)} for \eqn{\phi_{1} > \phi_{0}}{phi1 > phi0},
  or \code{northpole}\eqn{=c(0, 2*\pi - (\phi_{1}-\phi_{0})}{=c(0, 2\pi - (\phi_{1}-\phi_{0})} for \eqn{\phi_{1} < \phi_{0}}{phi1 < phi0}.
  More information about this rotation can be found in Section 3.2.2 of N. I. Fisher et al. (1993).

  Under this second rotation (i.e. where the colatitude of the point given
  in \code{northpole} has colatitude 0), the axis of rotation passes
  through the north and south poles as well as the origin, thus \eqn{x_1}{x1} is
  invariant under the rotation.  Note that a rotation where the colatitude
  of the point given in \code{northpole} is pi does not do the same, it
  sends the south pole to the north pole along the semicircle of longitude
  specified in \code{northpole}.  Thus for rotations under which the north
  and south poles are invariant, ensure that the colatitude of the point
  passed to \code{northpole} is 0. 

  If the points are to be rotated in a clockwise direction through a known
  angle about a known axis, then this can be done using \code{rot.sphere2},
  which applies Rodrigues' rotation formula (see, e.g. Koks, 2006).
  This can also be used if we know that the rotation will send a point \eqn{x}{x} to a point \eqn{y}{y} along
  the shortest great circle path.  In this case the axis of rotation, \eqn{v} is the
  unit length vector in the direction of the cross product of \eqn{x}{x} and \eqn{y}{y}, and \eqn{\theta=\pi}{theta=pi}.

  The function \code{rot.sphere2} can be used to rotate any point \eqn{x}{x} to the north pole \eqn{0_{sph}}{0sph}; simply let 
  \eqn{\theta}{theta} to be the colatitude of \eqn{x}{x} and \eqn{v}{v} to be the unit length vector in the direction of the cross product of \eqn{x}{x} and
  \eqn{0_{sph}}{0sph}.  To rotate a point to the south pole, set \eqn{v}{v} as for rotation to the north pole, and set
  \eqn{\theta}{theta} to be \eqn{\pi}{pi} + the colatitude of \eqn{x}{x}.

  To perform an inverse rotation using \code{rot.sphere2}, set \eqn{v}{v} to be as for the original rotation, and set
  \eqn{\theta}{theta} to be \eqn{2*\pi}{2pi} - the value of \eqn{\theta}{theta} used in the original rotation. 

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
  The output of \code{rot.sphere2} is in the same format as that which the argument code{points} is given in.
}
\references{
Fisher, N. I., Embleton, B. J. J. and Lewis, T. (1987) `Statistical Analysis of Spherical Data', 1st paperback edition.  Cambridge University Press.

Koks, D. (2006) `A Roundabout Route to Geometric Algebra', in \emph{Explorations in Mathematical Physics}, ch 4, pp. 144--157.  Springer Science + Business Media LLC.
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\examples{
rot.matrix(northpole=c(pi/3, pi/3), rad=1)

rot.sphere(points=c(pi/3, pi/3), northpole=c(pi/3,pi/3), rad=1, inverse=FALSE)
rot.sphere(points=c(pi/3, pi/3), northpole=c(pi/3,pi/3), rad=1, inverse=TRUE)

rot.sphere2(points=c(pi/3, pi/3), v=c(0,0), theta=pi)

}
\keyword{spatial}
\keyword{math}
