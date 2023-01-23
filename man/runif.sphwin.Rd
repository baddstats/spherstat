\name{runif.sphwin}
\alias{runif.sphwin}
\alias{rpoispp.sphwin}
\title{
  Generate the uniform distribution or a homogeneous Poisson point
  pattern on the sphere 
}
\description{
  Generates the uniform distribution or a random point pattern using the
  (homogeneous) Poisson process, on a sphere or a subset of a
  sphere. Includes CSR (complete spatial randomness).  
}
\usage{
runif.sphwin(n, win = sphwin(type = "sphere"), nsim=1, drop=TRUE,
 as.sp=TRUE, ndim="2")

rpoispp.sphwin(lambda, win = sphwin(type = "sphere"), lmax=NULL,
 ..., nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{n}{
    The number of points to be simulated.  An integer.
  }
  \item{lambda}{
    The intensity of the point process to be simulated.  A numeric of
    length 1 or an object of class \code{function}. 
  }
  \item{win}{
    The window, an object of type \code{sphwin}.
  }
  \item{lmax}{
    A numeric of length 1 giving the maximum value of \code{lambda}, if
    \code{lambda} is a function. 
  }
  \item{\dots}{
    Arguments passed to \code{lambda} if it is a function.
  }
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the result will be a point pattern, rather than a list containing a point pattern. 
  }
  \item{as.sp}{
    Logical.  If \code{TRUE}, returns an object of class as defined by
    \code{sp.dim}.  Otherwise, returns a matrix.  See Value. 
  }
  \item{ndim}{
    A string, taking value \code{"2"} or \code{"3"}.  Specifies whether the object should contain the locations of the points in spherical coordinates (\code{ndim="2"}) or Cartesian coordinates (\code{ndim="3"}).
  }

}
\details{
  \code{runif.sphwin} generates a realisation of \code{n} points of the
  uniform distribution inside the window \code{win} (if \code{win} is an
  object of class \code{sphwin}). 

  If \code{lambda} is a single number, then \code{rpoispp.sphwin} this
  algorithm generates a realisation of the uniform Poisson process (also
  known as Complete Spatial Randomness, CSR) inside the window
  \code{win} with intensity \code{lambda} (points per unit area).
  
  If \code{lambda} is a function, then \code{rpoispp.sphwin} uses the
  adaptation to point patterns on the sphere of Lewis and Shedler's (1979)
  rejection sampling algorithm.  According to this algorithm, a uniform
  Poisson process of intensity \code{lmax} is generated, and a point at
  location \code{s} is retained with probability \code{lambda(s)/lmax}.
}
\value{
  If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

  An item is determined by the values of \code{as.sp} and \code{ndim}:

  If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.

  If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{sp2} giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{sp3} giving the locations of the simulated points.
}
\references{
  Lewis, P.A.W. and Shedler, G.S. (1979)
  Simulation of non-homogeneous Poisson processes by thinning.
  \emph{Naval Logistics Quarterly} \bold{26}, 406--413.
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\note{
  The function \code{rpoispp.sphwin} is the analogue for point processes on
  the sphere of the function \code{\link[spatstat.random]{rpoispp}} in
  \pkg{spatstat}, which is the corresponding function for point
  processes in \eqn{R^2}{R^2}.  Hence elements of this help page have been taken
  from \code{\link[spatstat.random]{rpoispp}} with the permission of
  A. J. Baddeley.  This enables the relevant information on this help
  page to be consistent with that for \code{\link[spatstat.random]{rpoispp}}.
  It is hoped that this will minimise or remove any confusion for users
  of both \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link[spatstat.random]{rpoispp}},  \code{\link{rFisher}},
  \code{\link{rHardcore.sphwin}}, \code{\link{rMatClust.sphwin}},
  \code{\link{rMaternI.sphwin}}, \code{\link{rMaternII.sphwin}},
  \code{\link{rStrauss.sphwin}}, \code{\link{rThomas.sphwin}} 
}
\examples{
sph <- sphwin(type="sphere")
rus <- runif.sphwin(n=150, win=sph, as.sp=TRUE, ndim="2")
rus

rps <- rpoispp.sphwin(win=sph, lambda=150, as.sp=TRUE, ndim="2")
rps
}
\keyword{spatial}
\keyword{datagen}

