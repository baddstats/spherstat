\name{rotate.sphwin}
\alias{rotate.sphwin}
\alias{rotate.sp2}
\alias{rotate.sp3}
\title{
  Rotate a Spherical Window or Spherical Point Pattern
}
\description{
  These functions apply a rotation to a given spherical window
  or spherical point pattern.
}
\usage{
\method{rotate}{sphwin}(X, northpole, inverse = FALSE, \dots)

\method{rotate}{sp2}(X, northpole, inverse = FALSE, \dots)

\method{rotate}{sp3}(X, northpole, inverse = FALSE, \dots)
}
\arguments{
  \item{X}{
    Window or point pattern to be rotated: an object of class
    \code{"sphwin"}, \code{"sp2"} or \code{"sp3"}.
  }
  \item{northpole}{
    A point that will be rotated to the north pole (if
    \code{inverse=FALSE}, the default) or a point to which the north
    pole should be rotated (if \code{inverse=TRUE}).
    A numeric vector of length 2 containing
    polar coordinates (colatitude, longitude) in radians.
  }
  \item{inverse}{
    Logical value. If \code{FALSE} (the default), the rotation will be
    such that the point specified by \code{northpole} will move to the position
    of the north pole of the sphere. If \code{TRUE}, the rotation will
    be such that the north pole of the sphere will be rotated to the
    position given by \code{northpole}.
}
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions are methods for the generic
  \code{\link[spatstat]{rotate}}
  for the classes \code{"sphwin"}, \code{"sp2"} and \code{"sp3"}
  representing spherical windows and point patterns.

  A rotation on the sphere can be specified by stipulating that a point
  at position \code{A} should move to position \code{B}, where
  \code{A,B} are any given points on the sphere. The functions here
  assume that either \code{A} or \code{B} is the north pole.
  
  If \code{inverse=FALSE} (the default),
  the sphere is rotated so that the point specified by
  \code{northpole} moves to the north pole.
  If \code{inverse=TRUE}, the sphere is rotated so that the north pole
  is moved to the position specified by \code{northpole}.

  To apply a general rotation which moves point \code{A} to point
  \code{B}, use two rotations: \code{rotate(rotate(X, A), B, inverse=TRUE)}
}
\value{
  Another object of the same kind as \code{X} representing the result
  of rotation.
}
\author{
  Tom Lawrence and Adrian Baddeley
}
\examples{
  sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(pi/8,pi/6))  
  X <- runif.sphwin(30, sphcap)
  X
  Y <- rotate(X, c(0, pi/3))
  Y
}
\keyword{spatial}
\keyword{manip}
