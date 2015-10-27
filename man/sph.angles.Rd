\name{sph.angles}
\alias{sph.angles}
\title{
  Calculate angles in a spherical polygon
}
\description{
  Calculate the size of all angles in a spherical polygon
}
\usage{
  sph.angles(win)
}
\arguments{
  \item{win}{
    A window.  Must be an object of class \code{sphwin} with
    \code{type="polygon"}. 
  }
}
\details{
  This function uses the spherical geometry and the spherical cosine
  rule to calculate the sizes of all angles in a spherical polygon. 
}
\value{
  A numeric containing the sizes of all angles in radians, given in an
  anticlockwise direction.  The first angle in the numeric is at the
  second of the vertices in listed under \code{param} in the
  \code{sphwin} object. 
}
\author{
Tom Lawrence
}
\seealso{
\code{\link{sphcos}} (spherical cosine rule), \code{\link{sphsin}}
(spherical sine rule) 
}
\examples{
x <- matrix(nrow=4, ncol=2)
x[1,] <- c(0,0)
x[2,] <- c(pi/2, 0)
x[3,] <- c(pi/2, pi/2)
x[4,] <- c(0,0)
# sphpoly <- sphwin(type="polygon", param=x, ref=c(pi/4, pi/4))
sphpoly <- sphwin(type="polygon", param=x, ref=rep(0, 3))
sph.angles(win=sphpoly)
}
\keyword{spatial}
\keyword{math}

