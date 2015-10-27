\name{Kiso}
\alias{Kiso}
\title{
  Estimate of Spherical K Function
}
\description{
  Given a point pattern on a region of the sphere,
  compute an estimate of the \eqn{K}-function using isotropic edge correction.
}
\usage{
  Kiso(X, win, r, rad = win$rad, Dmat = pairdistsph(X),
       nrX = nrow(X), denom = (nrX * (nrX - 1))/area.sphwin(X))
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{K(r)} will be
    computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{win}{
    The window; an object of type \code{sphwin}.  Only required if
    \code{X} is an object of class \code{matrix}. 
  }
  \item{r}{
    Optional.  Numeric vector. The values of the argument \eqn{r} at
    which \code{K(r)} should be evaluated.  
  }
  \item{rad}{
    Sphere radius. A positive numeric value.
  }
  \item{Dmat}{
    Optional. Matrix of distances between all pairs of points in \code{X}.
    Distances are measured by great circle distance.
  }
  \item{nrX}{
    Optional. The number of points in \code{X}.
  }
  \item{denom}{
   Optional. The denominator for the \eqn{K} function.
  }
}
\details{
  To be written.
}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 
}
\author{
  Tom Lawrence.
}
\seealso{
  \code{\link{Fsphere}}
}
\examples{
   42
   # No examples for now - under development!
}
\keyword{spatial}
\keyword{nonparametric}
